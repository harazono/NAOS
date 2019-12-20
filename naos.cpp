#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdexcept>
#include "naos_library.h"
#include "cpas_debug.h"

#include "hts.h"
#include "sam.h"


using namespace std;


void printUsageAndExit(){
  cerr << "NAOS ref.fa align.sam" << endl;
  exit(2);
}
vector<alignment> alignments;
SequenceName      reference_name;


vector<allele> snv_candidates;

void parse_sam (
    const char* FASTAFileName,
    const char* SAMFileName,
    uint KmerSize,
    const bool outputInCSV,
    const string& binaryOutputFileName ///< empty() if --binary is not given
    )
{
  fprintf(stderr, "\n===Parameters===\n");
  fprintf(stderr, "Reference FASTA: %s\n", FASTAFileName);
  fprintf(stderr, "Input SAM file : %s\n", SAMFileName);
  fprintf(stderr, "K-mer size     : %d\n", KmerSize);
  fprintf(stderr, "================\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Loading from the FASTA file ...\n");
  const MultiFASTA multiFASTA = loadFromFASTA(FASTAFileName);
  fprintf(stderr, "Done                           \n");
  cerr << "reference record are as below\n---------" << endl;
  for(auto itr = multiFASTA.begin(); itr != multiFASTA.end(); itr++) {
    fprintf(stderr, "%s\n", itr->first.c_str());
  }
  cerr << "---------" << endl;
  FastTSVParse ftp(SAMFileName);
  if(!ftp) {
    fprintf(stderr, "ERROR: Cannot open SAM file '%s'\n", SAMFileName);
    exit(2);
  }

  SAMRecord record;
  size_t recordCount = 0;
  while(ftp.readNextLine()){
    const char* line = ftp.c_str();
    const bool isEmptyLine = line[0] == '\0';
    if(isEmptyLine) continue;
    const bool isCommentLine = line[0] == '@';
    if(isCommentLine) continue;
    if(!record.fill(ftp)) {
      cerr << "SAM Parsing ERROR at line " << ftp.getLineNumber();
      exit(2);
    }
    if(record.rname == "*") continue;
    if(record.seq == "*") continue;
    if(record.rname == record.qname) reference_name = record.rname;


    /*
     * use only primary alignment and supplementary alignment and reverse compliment of them.
     * sample has many supplimentary alignment.
     */
    if((record.flag & 2064) != record.flag) continue;

    if(!multiFASTA.count(record.rname.c_str())) {
      cerr << "SAM record says RNAME = '" << record.rname << "', but the reference genome does not have '" << record.rname << "'" << endl;
      continue;
      //exit(2);
    }

    // get aligned sequence at here
    const CIGAROPS cops     = parseCIGARString(record.cigar);

    BString refBS           = multiFASTA.at(record.rname);
    BString queryBS         = String2BString(record.seq);
    const int refStartPos   = record.pos;
    const int queryStartPos = 0; // generateAlignmentSequencesFromCIGARAndSeqs() will manege first Softclip / Hardclip
    BString ras, qas;
    generateAlignmentSequencesFromCIGARAndSeqs(refBS, queryBS, cops, refStartPos, queryStartPos, ras, qas);
    // if record flag has revcomp flag, modify aligned reference sequence and read sequence.

    /*  if((record.flag & 16) != 16){
        revCompBString(ras);
        revCompBString(qas);
        }
        */
    alignment al;
    al.ras           = ras;
    al.qas           = qas;
    al.ref_start_pos = refStartPos;
    al.readid        = recordCount;
    al.readname      = record.qname;
    alignments.push_back(al);
    //al.print_alignment();
    ++recordCount;
    if(recordCount % 100 == 0) {
      cerr << recordCount << " processed\r" << flush;
    }
  }
  cerr << recordCount << " processed\n";
  if(!binaryOutputFileName.empty()) {
    //outputAsBinaryTable(binaryOutputFileName);
  }
}


void pushback_snv_candidates(
    vector<alignment> const &alignments
    ){
  int recordCount = 0;
  alignment ref_al;
  fprintf(stdout, "enter to push back alignment\n");
  for(auto itr = alignments.begin(); itr != alignments.end(); itr++){
    if(itr->readname == reference_name) ref_al = *itr;
  }
  //ref_al.print_alignment();
  fprintf(stdout,"qas.size() = %d\n", BString2String(ref_al.ras).size());
  BString* ref_bstring = &ref_al.qas;
  string   ref_string  = BString2String(*ref_bstring);
  fprintf(stdout, "ref_string.size() = %d\n", ref_string.size());
  for(int ref_pos_itr = 0; ref_pos_itr < ref_string.size(); ref_pos_itr++){
    // check a position
    if(){
      continue;
    }

    fprintf(stdout, "ref_pos_itr = %d\n", ref_pos_itr);
    allele tmp_allele;
    if(ref_pos_itr < 10){
      tmp_allele.print_member();
    }
    tmp_allele.refallele = ref_string[ref_pos_itr];
    for(auto reads_itr = alignments.begin(); reads_itr != alignments.end(); reads_itr++) {
      fprintf(stdout,"read_itr->readid = %d\n", reads_itr->readid);
      tmp_allele.reads[reads_itr->readid] = reads_itr->ref_start_pos;
    }
    snv_candidates.push_back(tmp_allele);
    ++recordCount;
    if(recordCount % 100 == 0) {
      //cerr << recordCount << " processed\r" << flush;
    }
  }
  //cerr << recordCount << " processed\n";
}

void call_snv(
    allele al
    ){
}



int main(int argc, char *argv[]){
  GDB_On_SEGV g(argv[0]);
  struct option longopts[] = {
    { "csv"       , no_argument       , NULL , 'c' } ,
    // { "delete" , optional_argument , NULL , 'd' } ,
    { "kmer"      , required_argument , NULL , 'k' } ,
    { "binary"    , required_argument , NULL , 'b' } ,
    { 0           , 0                 , 0    , 0  }  ,
  };
  /// PARAMETERS ///

  int kmer_size = 1;
  bool output_in_csv = false;
  string binary_output_file_name;

  //////////////////

  int opt;
  int longindex;
  while ((opt = getopt_long(argc, argv, "k:", longopts, &longindex)) != -1) {
    switch (opt) {
      case 'b':
        binary_output_file_name = optarg;
        break;
      case 'c':
        output_in_csv = true;
        break;
      case 'k':
        kmer_size = atoi(optarg);
        if(kmer_size < 1 || 10 < kmer_size) {
          //error("kmer size must be 1-10");
        }
        break;
      default:
        MYASSERT_NEVERREACH();
    }
  }
  const int NUM_REQUIRED_ARGUMENTS = 2;
  if(optind + NUM_REQUIRED_ARGUMENTS != argc) {
    printUsageAndExit();
  }

  const char* fasta_file_name = argv[optind + 0];
  const char* sam_file_name   = argv[optind + 1];
  parse_sam(fasta_file_name, sam_file_name, kmer_size, output_in_csv, binary_output_file_name);
  pushback_snv_candidates(alignments);
  return 0;
}
