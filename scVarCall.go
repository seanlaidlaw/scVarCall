package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"strings"
	"sync"

	"github.com/spf13/viper"
)

// PIPELINE STEPS
// Sequencial steps
// 1. Split bam file to only MT chromosome
// 2. Deduplicate UMIs from bam file
// 3. Assess what barcodes are inside bam file

// Parallel steps
// 4. Split bam file by cell barcode
// 5. Samtools Quickcheck split bams
// 6. Sort bam file
// 7. Index bam file
// 8. Call all variants (not just second-max)

// fileExists checks if a file exists and is not a directory before we
// try using it to prevent further errors.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func writeCheckpoint(barcode_list []barcode_file, step int) {
	checkpoint_file := fmt.Sprintf("checkpoint_%d.json", step)
	rankingsJson, _ := json.MarshalIndent(barcode_list, "", "  ")
	err := ioutil.WriteFile(checkpoint_file, rankingsJson, 0644)
	if err != nil {
		panic(err)
	}
	log.Println(fmt.Sprintf("Checkpoint saved for step %d", step))
}

type barcode_file struct {
	Cell_barcode string
	Patient_ID   string
}

var barcode_list []barcode_file

var wg sync.WaitGroup

func main() {
	// we want to load a config file named "scVarCall.yaml" if it exists in WD or in ~/.config
	viper.SetConfigName("scVarCall")
	viper.SetConfigType("yaml")
	viper.AddConfigPath(".")              // look for config in the working directory first
	viper.AddConfigPath("$HOME/.config/") // if not found then look in .config folder

	viper.SetDefault("samtools_exec", "/software/sciops/pkgg/samtools/1.10.0/bin/samtools")
	viper.SetDefault("star_exec", "/nfs/users/nfs_r/rr11/Tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR")
	viper.SetDefault("star_genome_dir", "/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/75/")
	viper.SetDefault("genome_annot", "/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCH37d5/star/e75/ensembl.gtf")

	// read in config file if found, else use defaults
	if err := viper.ReadInConfig(); err != nil {
		log.Fatalln("Unable to read config file")
	}

	samtools_exec := viper.GetString("samtools_exec")

	var input string
	var output_dir string
	var threads string
	var current_step int

	// flags declaration using flag package
	flag.StringVar(&input, "i", "input", "input bam file produced by 10X CellRanger")
	flag.StringVar(&output_dir, "o", "output", "path to output directory")

	flag.Parse() // after declaring flags we need to call it
	if (strings.TrimSpace(input) == "input") || (strings.TrimSpace(output_dir) == "output_dir") {
		log.Fatalln("No input or output_dir argument was provided")
	}

	current_step = 0
	// if bam has been subset to MT already then skip to next step
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &barcode_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println("Subsetting bam file to MT only")
		err := os.Mkdir(output_dir, 0755)
		if err != nil {
			log.Fatal(err)
		}
		mt_subset_bam := output_dir + "/MT_subset.bam"
		output, err := exec.Command(
			"bsub",
			"-I",
			"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
			"-n", threads,
			samtools_exec, "view", "chrMT", "-b",
			"-@", threads, ">", mt_subset_bam).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		writeCheckpoint(barcode_list, current_step)
	}

}
