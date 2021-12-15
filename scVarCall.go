package main

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"reflect"
	"regexp"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/spf13/viper"
)

// PIPELINE STEPS
// Sequencial steps
// 0. Split bam file to only MT chromosome
// 0. Deduplicate UMIs from bam file
// 0. Assess what barcodes are inside bam file

// Parallel steps
// 0. Split bam file by cell barcode
// 0. Samtools Quickcheck split bams
// 0. Sort bam file
// 0. Index bam file
// 0. Call all variants (not just second-max)

func rmIfExists(file_path string) {
	if fileExists(file_path) {
		os.Remove(file_path)
	}
}

func chunkSlice(slice []barcode, chunkSize int) [][]barcode {
	var chunks [][]barcode
	for i := 0; i < len(slice); i += chunkSize {
		end := i + chunkSize

		// necessary check to avoid slicing beyond
		// slice capacity
		if end > len(slice) {
			end = len(slice)
		}

		chunks = append(chunks, slice[i:end])
	}

	return chunks
}

// fileExists checks if a file exists and is not a directory before we
// try using it to prevent further errors.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}
func stringInSlice(a string, list []string) bool {
	for _, b := range list {
		if b == a {
			return true
		}
	}
	return false
}

func writeCheckpoint(barcode_list []barcode, step int) {
	checkpoint_file := fmt.Sprintf("checkpoint_%d.json", step)
	checkpoint_file = output_dir + checkpoint_file
	rankingsJson, _ := json.MarshalIndent(barcode_list, "", "  ")
	err := ioutil.WriteFile(checkpoint_file, rankingsJson, 0644)
	if err != nil {
		panic(err)
	}
	log.Println(fmt.Sprintf("Checkpoint saved for step %d", step))
}

func bjobsIsCompleted(
	submitted_jobs_map map[string]string,
	attribute_name string,
	barcode_list *[]barcode,
	remove_list []string,
) {
	// while there are still jobs in submitted_jobs_map, iterate over reading through all job outputs and only leave when all have completed successfully
	for len(submitted_jobs_map) > 0 {
		for i := range *barcode_list {
			func_cram := &((*barcode_list)[i])
			bjobs_output_filename := submitted_jobs_map[func_cram.Name]

			dat, err := ioutil.ReadFile(bjobs_output_filename)
			if err == nil {
				// when job has finished (either successfully or with exit code, remove from the waiting list 'submitted_jobs_map'
				if strings.Contains(string(dat), "Terminated at") {
					delete(submitted_jobs_map, func_cram.Name)

					// if job has finished and successfully completed then set the specified attribute_name to true
					if strings.Contains(string(dat), "Successfully completed.") {
						reflect.ValueOf(func_cram).Elem().FieldByName(attribute_name).SetBool(true)

						if len(remove_list) > 0 {
							for _, file := range remove_list {
								filename_path := reflect.ValueOf(func_cram).Elem().FieldByName(file).Interface().(string)
								err := os.Remove(filename_path)
								if err != nil {
									log.Fatal(err)
								}
							}
						}

					} else {
						log.Println(fmt.Sprintf("Error with bsub job: %s", bjobs_output_filename))
					}
				}
			}
		}
		// sleep for 5 seconds after going through every job's output before retrying
		time.Sleep(5 * time.Second)
	}
}

//func quickcheck_alignments(barcode_list []barcode, i int, samtools_exec string) {
//	cram := &barcode_list[i]
//
//	bam_filename := cram.Realigned_bam_path
//	output, err := exec.Command(samtools_exec, "quickcheck", bam_filename).CombinedOutput()
//
//	if err != nil {
//		// Display everything we got if error.
//		log.Println("Error when running command.  Output:")
//		log.Println(string(output))
//		log.Printf("Got command status: %s\n", err.Error())
//		cram.Realigned_quickcheck_success = false
//	} else {
//		cram.Realigned_quickcheck_success = true
//	}
//	wg.Done()
//}

func indexBam(bam_filename string) {
	output, err := exec.Command(samtools_exec, "index", bam_filename).CombinedOutput()

	if err != nil {
		// Display everything we got if error.
		log.Println("Error when running command.  Output:")
		log.Println(string(output))
		log.Printf("Got command status: %s\n", err.Error())
	}
}

type barcode struct {
	Name                                     string
	Output_dir                               string
	Masterbam_original                       string
	Masterbam_original_quickcheck_success    bool
	Masterbam_MT_subset                      string
	Masterbam_MT_subset_quickcheck_success   bool
	Masterbam_MT_subset_index_success        bool
	Masterbam_QC_subset                      string
	Masterbam_QC_subset_quickcheck_success   bool
	Masterbam_QC_subset_index_success        bool
	Masterbam_UMI_deduped                    string
	Masterbam_UMI_deduped_success            bool
	Masterbam_UMI_deduped_quickcheck_success bool
	Masterbam_UMI_deduped_index_success      bool
	Splitbam_jobout                          string
	Splitbam_joberr                          string
	Splitbam_bamout                          string
	Splitbam_bamindex                        string
	Splitbam_barcodefile                     string
	Splitbam_successful                      bool
	Splitbam_indexed                         bool
	Rvarcall_jobout                          string
	Rvarcall_joberr                          string
	Rvarcall_dir_out                         string
	Rvarcall_call_out                        string
	Rvarcall_cov_out                         string
	Rvarcall_command_successful              bool
	Rdsmerge_rds_calls                       string
	Rdsmerge_rds_coverage                    string
	Rdsmerge_success                         bool
}

var barcode_list []barcode
var chunk_cell_map map[string]string

var output_dir string
var samtools_exec string
var Rscript_exec string

var wg sync.WaitGroup

func main() {
	// we want to load a config file named "scVarCall.yaml" if it exists in WD or in ~/.config
	viper.SetConfigName("scVarCall")
	viper.SetConfigType("yaml")
	viper.AddConfigPath(".")              // look for config in the working directory first
	viper.AddConfigPath("$HOME/.config/") // if not found then look in .config folder

	viper.SetDefault("samtools_exec", "/software/sciops/pkgg/samtools/1.10.0/bin/samtools")
	viper.SetDefault("Rscript_exec", "/software/R-4.1.0/bin/Rscript")
	viper.SetDefault("umitools_exec", "/software/teamtrynka/conda/trynka-base/bin/umi_tools")
	viper.SetDefault(
		"star_exec",
		"/nfs/users/nfs_r/rr11/Tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR",
	)
	viper.SetDefault(
		"star_genome_dir",
		"/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/75/",
	)
	viper.SetDefault(
		"genome_annot",
		"/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCH37d5/star/e75/ensembl.gtf",
	)

	// read in config file if found, else use defaults
	//if err := viper.ReadInConfig(); err != nil {
	//	log.Fatalln("Unable to read config file")
	//}

	samtools_exec = viper.GetString("samtools_exec")
	Rscript_exec = viper.GetString("Rscript_exec")
	umitools_exec := viper.GetString("umitools_exec")
	//star_exec := viper.GetString("star_exec")
	//star_genome_dir := viper.GetString("star_genome_dir")
	//featurecounts_exec := viper.GetString("featurecounts_exec")
	//genome_annot := viper.GetString("genome_annot")

	var input string
	var barcodes_qc string
	var threads string
	var current_step int

	// flags declaration using flag package
	flag.StringVar(&input, "i", "input", "input bam file produced by 10X CellRanger")
	flag.StringVar(&output_dir, "o", "output", "path to output directory")
	flag.StringVar(&barcodes_qc, "b", "barcodes", "list of QC passed barcodes")

	flag.Parse() // after declaring flags we need to call it
	if (strings.TrimSpace(input) == "input") || (strings.TrimSpace(output_dir) == "output_dir") {
		log.Fatalln("No input or output_dir argument was provided")
	}
	// make sure output dir ends in slash so paths work correctly when appending filenames
	output_dir = output_dir + "/"

	mt_subset_bam := output_dir + "/MT_subset.bam"
	threads = "4"

	current_step = 1
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Defining input and output paths for master bam")
		err := os.Mkdir(output_dir, 0755)
		if err != nil {
			log.Fatal(err)
		}

		// define a "master" barcode for the bam file that is to be split
		// this isnt a barcode per se, but the origin of the barcodes
		var master_barcode barcode
		master_barcode.Name = "MASTER"
		master_barcode.Output_dir = output_dir
		master_barcode.Masterbam_original = input

		barcode_list = append(barcode_list, master_barcode)
		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 2
	// if bam has been subset to MT already then skip to next step
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		master_barcode := &barcode_list[0]

		log.Println("Quickchecking input bam file")

		output, err := exec.Command(samtools_exec, "quickcheck", master_barcode.Masterbam_original).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			master_barcode.Masterbam_original_quickcheck_success = false
		} else {
			master_barcode.Masterbam_original_quickcheck_success = true
		}

		log.Println("Subsetting bam file to MT only")

		output, err = exec.Command(
			"bsub",
			"-I",
			"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
			"-n", threads,
			samtools_exec, "view",
			input, "MT",
			"-b", "-@", threads,
			">", mt_subset_bam).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		master_barcode.Masterbam_MT_subset = mt_subset_bam

		log.Println("Quickchecking subset bam file")

		output, err = exec.Command(samtools_exec, "quickcheck", master_barcode.Masterbam_MT_subset).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			master_barcode.Masterbam_MT_subset_quickcheck_success = false
		} else {
			master_barcode.Masterbam_MT_subset_quickcheck_success = true
		}

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 3
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Indexing newly created MT subset bam")

		indexBam((&barcode_list[0]).Masterbam_MT_subset)
		(&barcode_list[0]).Masterbam_MT_subset_index_success = true

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 4
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Subsetting to QC passed barcodes")
		(&barcode_list[0]).Masterbam_QC_subset = output_dir + "/MT_subset_QC_filtered.bam"

		output, err := exec.Command(
			"bsub",
			"-I",
			"-R'select[mem>5000] rusage[mem=5000]'", "-M5000",
			"-n", "12",
			"subset-bam", "--cores", "12",
			"--bam", (&barcode_list[0]).Masterbam_MT_subset,
			"--cell-barcodes", barcodes_qc,
			"--out-bam", (&barcode_list[0]).Masterbam_QC_subset).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		output, err = exec.Command(samtools_exec, "quickcheck", (&barcode_list[0]).Masterbam_QC_subset).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			(&barcode_list[0]).Masterbam_QC_subset_quickcheck_success = false
		} else {
			(&barcode_list[0]).Masterbam_QC_subset_quickcheck_success = true
		}

		indexBam((&barcode_list[0]).Masterbam_QC_subset)
		(&barcode_list[0]).Masterbam_QC_subset_index_success = true

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 5
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Deduplicating UMIs")
		deduped_bam := output_dir + "/MT_subset_umi_deduped.bam"
		(&barcode_list[0]).Masterbam_UMI_deduped = deduped_bam

		output, err := exec.Command(
			"bsub",
			"-I",
			"-R'select[mem>80000] rusage[mem=80000]'", "-M80000",
			umitools_exec, "dedup",
			"--paired",
			"--chrom", "MT",
			"--extract-umi-method", "tag",
			"--umi-tag", "UB",
			"--per-cell", "--cell-tag", "CB",
			"-I", (&barcode_list[0]).Masterbam_QC_subset,
			"-S", (&barcode_list[0]).Masterbam_UMI_deduped).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		(&barcode_list[0]).Masterbam_UMI_deduped_success = true

		// quickcheck produced bam
		output, err = exec.Command(samtools_exec, "quickcheck", (&barcode_list[0]).Masterbam_UMI_deduped).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			(&barcode_list[0]).Masterbam_UMI_deduped_quickcheck_success = false
		} else {
			(&barcode_list[0]).Masterbam_UMI_deduped_quickcheck_success = true
		}

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 6
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Indexing newly created deduped bam")

		indexBam((&barcode_list[0]).Masterbam_UMI_deduped)
		(&barcode_list[0]).Masterbam_UMI_deduped_index_success = true

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 7
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))
		log.Println("Reading barcodes in deduped and subset bam input")
		output, err := exec.Command(
			"bsub",
			"-I",
			"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
			"-n", threads,
			samtools_exec, "view", (&barcode_list[0]).Masterbam_UMI_deduped,
			"|", "grep", "-oE", "CB:Z:[acgtnACGTN-]+[1-9]",
			"|", "sort", "-u", "--parallel", threads,
			"|", "gzip", ">", output_dir+"unique_barcodes.tsv.gz").CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		writeCheckpoint(barcode_list, current_step)
	}

	current_step = 8
	if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		//log.Println("Parsing unique barcodes file into memory")
		barcodefile, err := os.Open(output_dir + "unique_barcodes.tsv.gz")
		if err != nil {
			log.Fatal(err)
		}
		defer barcodefile.Close()

		gr, err := gzip.NewReader(barcodefile)
		if err != nil {
			log.Fatal(err)
		}
		defer gr.Close()

		barcodefile_scanner := bufio.NewScanner(gr)
		for barcodefile_scanner.Scan() {
			// remove 'CB:Z:' prefix to get just the cell barcode
			leading_regex := regexp.MustCompile(`^CB:Z:`)
			barcode_str := leading_regex.ReplaceAllString(barcodefile_scanner.Text(), "")

			if len(barcode_str) != 18 {
				log.Println("Problem with barcode parsed from: " + barcodefile_scanner.Text())
				log.Fatalln(fmt.Sprintf("Barcode '%s' is not of expected length, should be 18bp (with -1 ending) but is %d", barcode_str, len(barcode_str)))
			}

			new_barcode := barcode{Name: barcode_str}

			barcode_list = append(barcode_list, new_barcode)
		}

		// chunk the list of barcodes into groups of 500
		chunked_barcode_list := chunkSlice(barcode_list[1:], 500)

		log.Println("Spliting and calling variants on chunks of 500 barcodes")
		for chunk_i, chunk := range chunked_barcode_list {

			chunk_cell_map = make(map[string]string)
			chunk_output := output_dir + "/chunk_" + strconv.Itoa(chunk_i) + "/"
			err := os.Mkdir(chunk_output, 0755)
			if err != nil {
				log.Fatal(err)
			}

			for i := range chunk {
				cell := &chunk[i]

				cell.Splitbam_jobout = chunk_output + "cellsplit_" + cell.Name + ".o"
				cell.Splitbam_joberr = chunk_output + "cellsplit_" + cell.Name + ".e"
				cell.Splitbam_bamout = chunk_output + "cell_" + cell.Name + ".bam"
				cell.Splitbam_bamindex = cell.Splitbam_bamout + ".bai"
				cell.Splitbam_barcodefile = chunk_output + "barcode_" + cell.Name + ".txt"

				// write a file with barcode inside for splitbam
				job_barcode_out, err := os.Create(cell.Splitbam_barcodefile)
				if err != nil {
					log.Fatal(err)
				}
				defer job_barcode_out.Close()
				_, err = job_barcode_out.WriteString(cell.Name + "\n")

				// add this barcode to list of current jobs
				chunk_cell_map[cell.Name] = cell.Splitbam_jobout

				// split bam to current barcode
				output, err := exec.Command(
					"bsub",
					"-o", cell.Splitbam_jobout,
					"-e", cell.Splitbam_joberr,
					"-R'select[mem>5000] rusage[mem=5000]'", "-M5000",
					"-n", "12",
					"subset-bam", "--cores", "12",
					"--bam", (&barcode_list[0]).Masterbam_UMI_deduped,
					"--cell-barcodes", cell.Splitbam_barcodefile,
					"--out-bam", cell.Splitbam_bamout).CombinedOutput()

				if err != nil {
					// Display everything we got if error.
					log.Println("Error when running command.  Output:")
					log.Println(string(output))
					log.Printf("Got command status: %s\n", err.Error())
					return
				}
			}

			// wait for that chunks splits to finish
			bjobsIsCompleted(chunk_cell_map, "Splitbam_successful", &barcode_list, []string{"Splitbam_jobout", "Splitbam_joberr", "Splitbam_barcodefile"})

			// index the split bam files
			for i := range chunk {
				cell := &chunk[i]
				if cell.Splitbam_successful {
					indexBam(cell.Splitbam_bamout)
					cell.Splitbam_indexed = true
				}
			}

			// run variant calling on the bam files
			for i := range chunk {
				cell := &chunk[i]
				if cell.Splitbam_indexed {

					cell.Rvarcall_jobout = chunk_output + "Rvarcall_" + cell.Name + ".o"
					cell.Rvarcall_joberr = chunk_output + "Rvarcall_" + cell.Name + ".e"
					cell.Rvarcall_dir_out = chunk_output

					// add this barcode to list of current jobs
					chunk_cell_map[cell.Name] = cell.Rvarcall_jobout

					// call variants on bam
					output, err := exec.Command(
						"bsub",
						"-o", cell.Rvarcall_jobout,
						"-e", cell.Rvarcall_joberr,
						"-R'select[mem>5000] rusage[mem=5000]'", "-M5000",
						"-n", "1",
						Rscript_exec, "callVars.R",
						cell.Splitbam_bamout,
						cell.Rvarcall_dir_out).CombinedOutput()

					if err != nil {
						// Display everything we got if error.
						log.Println("Error when running command.  Output:")
						log.Println(string(output))
						log.Printf("Got command status: %s\n", err.Error())
						return
					}

					// set expected output merged rds filenames to barcode object
					barcode_trimmed := strings.ReplaceAll(cell.Name, "-1", "")
					cell.Rvarcall_call_out = cell.Rvarcall_dir_out + "cell_" + barcode_trimmed + ".calls.rds"
					cell.Rvarcall_cov_out = cell.Rvarcall_dir_out + "cell_" + barcode_trimmed + ".coverage.rds"

				}
			}

			// wait for variant calls to finish
			bjobsIsCompleted(chunk_cell_map, "Rvarcall_command_successful", &barcode_list, []string{"Rvarcall_jobout", "Rvarcall_joberr", "Splitbam_bamout", "Splitbam_bamindex"})

			// merge completed rds together into one rds for whole chunk
			output, err := exec.Command(
				"bsub",
				"-I",
				"-R'select[mem>16000] rusage[mem=16000]'", "-M16000",
				Rscript_exec, "mergeVarcallRds.R",
				chunk_output, strconv.Itoa(chunk_i)).CombinedOutput()

			if err != nil {
				// Display everything we got if error.
				log.Println("Error when running command.  Output:")
				log.Println(string(output))
				log.Printf("Got command status: %s\n", err.Error())
				return
			}

			for i := range chunk {
				cell := &chunk[i]
				if cell.Rvarcall_command_successful {

					cell.Rdsmerge_rds_coverage = chunk_output + "chunk_" + strconv.Itoa(chunk_i) + ".coverage.rds"
					cell.Rdsmerge_rds_calls = chunk_output + "chunk_" + strconv.Itoa(chunk_i) + ".calls.rds"

					cell.Rdsmerge_success = true

					if fileExists(cell.Rdsmerge_rds_coverage) {
						rmIfExists(cell.Rvarcall_cov_out)
					}

					if fileExists(cell.Rdsmerge_rds_calls) {
						rmIfExists(cell.Rvarcall_call_out)
					}
				}

			}
		}
		writeCheckpoint(barcode_list, current_step)

	}

	//current_step = 8
	//if fileExists(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(output_dir + fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {
	//log.Println(fmt.Sprintf("Starting step %d", current_step))

	//log.Println("Merging the Rds tables by chunk")

	//chunk_rds_map := make(map[string]string)

	//var chunk_list []string
	//for _, cell := range barcode_list[1:] {
	//if !stringInSlice(cell.Rvarcall_out, chunk_list) {
	//chunk_list = append(chunk_list, cell.Rvarcall_out)
	//}
	//}

	//for chunk_i, chunk := range chunk_list {
	//chunk_output := output_dir + "/chunk_" + strconv.Itoa(chunk_i) + "/"

	//rdsmerge_out := chunk_output + "rdsmerge.out"
	//rdsmerge_err := chunk_output + "rdsmerge.err"

	//output, err := exec.Command(
	//"bsub",
	//"-o", rdsmerge_out,
	//"-e", rdsmerge_err,
	//"-R'select[mem>16000] rusage[mem=16000]'", "-M16000",
	//Rscript_exec, "mergeVarcallRds.R",
	//chunk, strconv.Itoa(chunk_i)).CombinedOutput()

	//if err != nil {
	//// Display everything we got if error.
	//log.Println("Error when running command.  Output:")
	//log.Println(string(output))
	//log.Printf("Got command status: %s\n", err.Error())
	//return
	//}

	//for i := range barcode_list[1:] {
	//cell := &barcode_list[i]
	//if cell.Rvarcall_out == chunk {
	//chunk_rds_map[cell.Name] = rdsmerge_out
	//cell.Rdsmerge_out = rdsmerge_out

	//cell.Rdsmerge_err = rdsmerge_err
	//cell.Rdsmerge_rds_coverage = chunk_output + "chunk_" + strconv.Itoa(chunk_i) + ".coverage.rds"
	//cell.Rdsmerge_rds_calls = chunk_output + "chunk_" + strconv.Itoa(chunk_i) + ".calls.rds"
	//}
	//}

	//}
	//bjobsIsCompleted(chunk_rds_map, "Rdsmerge_success", &barcode_list, []string{})
	//writeCheckpoint(barcode_list, current_step)
	//}

	//// if symlinks already created then load session information and continue
	//// if fastq have already been split then load checkpoint instead of rerunning
	//if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {

	//log.Println("Symlinking fastq into different folders based on library_type")
	//for i := range barcode_list {
	//cram := &barcode_list[i]
	//if cram.Fastq_extracted_success && cram.Library_type != "" {
	//lib_type_dir := strings.ReplaceAll(cram.Library_type, " ", "_")
	//lib_type_dir = "C_Split_by_Library_Type/" + lib_type_dir
	//_ = os.MkdirAll(lib_type_dir, 0755)
	//os.Symlink("../../"+cram.Fastq_1_path, lib_type_dir+"/"+cram.Sample_name+".1.fq.gz")
	//os.Symlink("../../"+cram.Fastq_2_path, lib_type_dir+"/"+cram.Sample_name+".2.fq.gz")
	//cram.Symlinked_fq_1 = lib_type_dir + "/" + cram.Sample_name + ".1.fq.gz"
	//cram.Symlinked_fq_2 = lib_type_dir + "/" + cram.Sample_name + ".2.fq.gz"
	//}
	//}

	//writeCheckpoint(barcode_list, current_step)
	//}

	//current_step = 7
	//// if alignments already performed then load session information and continue
	//if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {

	//_ = os.MkdirAll("D_realignments/", 0755)
	//log.Println("Running alignments between extracted fastq and specified reference")
	//realignment_map := make(map[string]string)
	//for i := range barcode_list {
	//cram := &barcode_list[i]
	//if cram.Symlinked_fq_1 != "" && cram.Symlinked_fq_2 != "" {

	//out_folder := "D_realignments/" + strings.ReplaceAll(cram.Library_type, " ", "_") + "/"
	//_ = os.MkdirAll(out_folder, 0755)

	//bam_output := out_folder + cram.Sample_name + ".bam"
	//job_out := out_folder + "/D_realignement_RNA_" + cram.Sample_name + ".o"
	//job_err := out_folder + "/D_realignement_RNA_" + cram.Sample_name + ".e"

	//if stringInSlice(cram.Library_type, star_align_libraries) {
	//output, err := exec.Command(
	//"bsub",
	//"-o", job_out,
	//"-e", job_err,
	//"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
	//"-n", "10",
	//star_exec, "--runThreadN", "10",
	//"--outSAMattributes", "NH", "HI", "NM", "MD",
	//"--limitBAMsortRAM", "31532137230",
	//"--outSAMtype", "BAM", "SortedByCoordinate",
	//"--genomeDir", star_genome_dir,
	//"--readFilesCommand", "zcat",
	//"--outFileNamePrefix", out_folder+"/"+cram.Filename,
	//"--readFilesIn", cram.Symlinked_fq_1, cram.Symlinked_fq_2,
	//"--outStd", "BAM_SortedByCoordinate",
	//"|", samtools_exec, "sort", "-@3", "-l7", "-o", bam_output).CombinedOutput()

	//if err != nil {
	//// Display everything we got if error.
	//log.Println("Error when running command.  Output:")
	//log.Println(string(output))
	//log.Printf("Got command status: %s\n", err.Error())
	//return
	//}

	//cram.Realigned_bam_path = bam_output
	//realignment_map[cram.Filename] = job_out

	//} else if stringInSlice(cram.Library_type, bwa_align_libraries) {
	//output, err := exec.Command(
	//"bsub",
	//"-o", job_out,
	//"-e", job_err,
	//"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
	//"-n", "10",
	//bwa_exec, "mem", "-t", "10",
	//bwa_genome_ref,
	//cram.Symlinked_fq_1,
	//cram.Symlinked_fq_2,
	//"|", samtools_exec, "sort", "-@3", "-l7", "-o", bam_output).CombinedOutput()

	//if err != nil {
	//// Display everything we got if error.
	//log.Println("Error when running command.  Output:")
	//log.Println(string(output))
	//log.Printf("Got command status: %s\n", err.Error())
	//return
	//}

	//cram.Realigned_bam_path = bam_output
	//realignment_map[cram.Filename] = job_out
	//}
	//}
	//}

	//// check on alignment jobs until they have finished
	//bjobsIsCompleted(realignment_map, "Realigned_succesful", &barcode_list)

	//writeCheckpoint(barcode_list, current_step)
	//}

	//current_step = 8
	//// if quickcheck has already been performed then load session information and continue
	//if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {
	//log.Println("Running samtools quickcheck on completed bams")

	//for i := range barcode_list {
	//cram := &barcode_list[i]
	//if cram.Realigned_succesful {
	//wg.Add(1)
	//go quickcheck_alignments(barcode_list, i, samtools_exec)
	//}
	//}
	//wg.Wait() // wait until all quickcheck processes have finished

	//writeCheckpoint(barcode_list, current_step)
	//}

	//current_step = 9
	//// if quickcheck has already been performed then load session information and continue
	//if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {
	//log.Println("Indexing quickchecked bams")

	//for i := range barcode_list {
	//cram := &barcode_list[i]
	//if cram.Realigned_quickcheck_success {
	//wg.Add(1)
	//indexBam(barcode_list, i, samtools_exec)
	//}
	//}
	//wg.Wait() // wait until all quickcheck processes have finished

	//writeCheckpoint(barcode_list, current_step)
	//}

	//current_step = 10
	//// if counts matrix has already been built then load session info
	//if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
	//jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
	//byteValue, _ := ioutil.ReadAll(jsonFile)
	//err = json.Unmarshal([]byte(byteValue), &barcode_list)
	//if err != nil {
	//panic(err)
	//}

	//log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	//} else {
	//log.Println("Running featurecounts on completed RNA bams")
	//var rna_bams_featurecounts_input []string

	//for i := range barcode_list {
	//cram := &barcode_list[i]
	//// if quickcheck worked then add its realigned and sorted bam path to list of bams to include in counts matrix
	//if cram.Realigned_quickcheck_success {
	//if stringInSlice(cram.Library_type, star_align_libraries) {
	//rna_bams_featurecounts_input = append(rna_bams_featurecounts_input, cram.Realigned_bam_path)

	//}
	//}
	//}

	//if len(rna_bams_featurecounts_input) < 1 {
	//log.Fatalln("Less than  1 bams in RNA category, not enough for featurecounts, aborting.")
	//}

	//err := os.Mkdir("E_Counts_matrix_RNA", 0755)
	//if err != nil {
	//log.Fatal(err)
	//}

	//matrix_out := "E_Counts_matrix_RNA/featurecounts_matrix.tsv"
	//job_out := "E_Counts_matrix_RNA/featurecounts_run.o"
	//job_err := "E_Counts_matrix_RNA/featurecounts_run.e"

	//featureCountsCmd := []string{
	//"-o", job_out,
	//"-e", job_err,
	//"-R'select[mem>20000] rusage[mem=20000]'", "-M20000",
	//"-n", "14",
	//featurecounts_exec,
	//"-Q", "30",
	//"-p",
	//"-t", "exon",
	//"-g", "gene_name",
	//"-F", "GTF",
	//"-a", genome_annot,
	//"-o", matrix_out}

	//// append bam paths to end of command options, as this is what featureCounts expects
	//featureCountsCmd = append(featureCountsCmd, rna_bams_featurecounts_input...)

	//output, err := exec.Command("bsub", featureCountsCmd...).CombinedOutput()

	//if err != nil {
	//// Display everything we got if error.
	//log.Println("Error when running command.  Output:")
	//log.Println(string(output))
	//log.Printf("Got command status: %s\n", err.Error())
	//return
	//}

	//// wait for featurecounts job to finish
	//for {
	//dat, err := ioutil.ReadFile(job_out)
	//if err == nil {
	//if strings.Contains(string(dat), "Terminated at") {
	//break
	//}
	//}
	//time.Sleep(5 * time.Second)
	//}

	//// if featurecounts exited successfuly write new checkpoint file
	//// this doesn't have any new information but its presence will indicate not to repeat the featurecounts step
	//dat, err := ioutil.ReadFile(job_out)
	//if err == nil {
	//if strings.Contains(string(dat), "Successfully completed.") {
	//writeCheckpoint(barcode_list, current_step)
	//}
	//}
	//}
}
