#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define arraylen(x) (sizeof (x) / sizeof (*x))
#define newline (printf("\n"))
#define douline (printf("\n\n"))
#define MASTER_NODE 0
#define SLAVE_NODE 1
#define add_Str_Terminator '\0'
#define SHELLSCRIPT "sudo chmod 755 movefile.sh && ./movefile.sh 2> /dev/null"

/*
 define constant variables for standard S. Pyogenes Cas9.
 For a more flexible program these could be set as function parameters.
*/
typedef enum {
	GUIDE_LEN = 20, 
	PAM_LEN = 3, // protospacer adjacent motif (2 - 6)
	CUT_DIFF = 3
} parameters;

/*
 Helper function that returns the compliment of a DNA sequence.
 */
char dna_compliment(char nucleotide) {
    
    char compliment;

    if (nucleotide == 'A') {
        compliment = 'T';
    } else if (nucleotide == 'T') {
        compliment = 'A';
    } else if (nucleotide == 'C') {
        compliment = 'G';
    } else if (nucleotide == 'G') {
        compliment = 'C';
    } else {
        compliment = '\0';
    }
    
    return compliment;
}

/*
 Helper function that makes sure the program does not go out of bounds.
 Trying to find target sequences at the start or end of the sequence.
 */
int enough_seq_context(int pam_start, int seq_len, char* strand) {

    if (strncmp(strand, "FWD", 3) == 0 && (pam_start > GUIDE_LEN && pam_start + PAM_LEN) <= seq_len) {
        return 1;
    }
    if (strncmp(strand, "REV", 3) == 0 && (pam_start + PAM_LEN + GUIDE_LEN) < seq_len) {
        return 1;
    }
    
    return 0;
}

/*
 Helper function that checks if sequence matches PAM.
 */
int pam_matches(char* pam_sequence) {
    
    static const char PAM[PAM_LEN] = "NGG"; // N is any nucleotide (the ‘letters’ that make up DNA)
    
    if (strlen(pam_sequence) != PAM_LEN) { return 0; }

    for (int i = 0; i < PAM_LEN; i++) {
        if (PAM[i] == 'N') {
            //pass (skip)
        }
        else if (*(pam_sequence + i) != PAM[i]) {
            return 0;
        }
    }
    
    return 1;
}

/*
 Helper function that read txt file as program input.
 */
char *readFile(char *fileName) {
    
    FILE *fp = fopen(fileName, "r");
    char *dna_file;
    size_t n = 0;
    int C;
    
    if (fp == NULL) {
        return NULL; //could not open file
    }
    
    fseek(fp, 0, SEEK_END);
    long f_size = ftell(fp); //get file size
    fseek(fp, 0, SEEK_SET);

    dna_file = malloc(f_size * sizeof(char)); //declare malloc size of read file
    
    while ((C = fgetc(fp)) != EOF) {
        dna_file[n++] = (char)C; //append charactor by charactor from read file into dna_file
    }
    
    dna_file[n] = '\0'; //insert \0 as the last charactor to notify end of string

    fclose(fp);
    
    return dna_file;
}

/*
 Helper function that returns the reverse compliment of a DNA sequence.
 Useful for returning results from the non reference (bottom) strand
 while still keeping the 5' to 3' convention.
 Implement XOR to swap element from edge to center.
 */
char *reverse(char* strand) {

    char *p1, *p2;

    if (! strand || ! *strand) { return strand; }

    for (p1 = strand, p2 = strand + strlen(strand) - 1; p2 > p1; p1++, p2--) {
        *p1 ^= *p2;
        *p2 ^= *p1;
        *p1 ^= *p2;
    }

    return strand;
}

int main(int argc, char *argv[])
{   
    char* input_strand = malloc(sizeof(char));
    input_strand = readFile("strand.txt");
    
    // Sequence length
    const int seq_len = (int) strlen(input_strand);

    int max_element = seq_len / PAM_LEN;

    /* MPI declare variables */
    int rank, size, error_code;
    double str_time, end_time, execTime, elapsedTime;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Quit MPI program if processors != 2
    if (size != 2) { 
        printf("Quitting. Number of processors must be 2.\n");
        MPI_Abort(MPI_COMM_WORLD, error_code);
        exit(0);
    }

    // Get processors name
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;

    MPI_Barrier(MPI_COMM_WORLD);

    str_time = MPI_Wtime(); // start timer

    // input_strand --> sequence (pointer)
    char *sequence;
    MPI_Alloc_mem((seq_len * (sizeof *sequence) + sizeof(char)), MPI_INFO_NULL, &sequence);
    memcpy(sequence, input_strand, seq_len + 1); // included '\0'
    
    int *cut_pos;
    MPI_Alloc_mem(max_element * (sizeof *cut_pos), MPI_INFO_NULL, &cut_pos);
    char (*pam_seq)[PAM_LEN + 1];
    MPI_Alloc_mem(max_element * (sizeof *pam_seq), MPI_INFO_NULL, &pam_seq);
    char (*target_seq)[GUIDE_LEN + 1];
    MPI_Alloc_mem(max_element * (sizeof *target_seq), MPI_INFO_NULL, &target_seq);
    
    char pam_seq_temp[PAM_LEN + 1];
    int local_idx = 0;

    if (rank == MASTER_NODE)
    {   
        MPI_Get_processor_name(processor_name, &name_len);

        for (int i = GUIDE_LEN + CUT_DIFF; i < seq_len; i++) {

            for (int j = i, idx = 0; j < i + PAM_LEN; j++, idx++) {
                pam_seq_temp[idx] = *(sequence + j);
            }
            pam_seq_temp[PAM_LEN] = add_Str_Terminator; // add '\0' at the end

            // find targets on reference (top) strand
            // check for out of bounds and correct PAM
            if (enough_seq_context(i, seq_len, "FWD") && pam_matches(pam_seq_temp)) {
                // We take the reverse compliment of the pam_seq and target_seq to
                // keep the 5` to 3` convention

                cut_pos[local_idx] = i - CUT_DIFF;
                    
                for (int j = i, pam_sub = 0; j < i + PAM_LEN; j++, pam_sub++) {
                    pam_seq[local_idx][pam_sub] = *(sequence + j);
                }
                pam_seq[local_idx][PAM_LEN] = add_Str_Terminator;
                    
                for (int j = i - GUIDE_LEN, tar_sub = 0; j < i; j++, tar_sub++) {
                    target_seq[local_idx][tar_sub] = *(sequence + j);
                }
                target_seq[local_idx][GUIDE_LEN] = add_Str_Terminator;

                local_idx++;
            }
        }

        // write output to file forward.txt
        FILE *fpr = fopen("forward.txt", "wb");

        if (fpr != NULL) {

            fputs("Programming the genome with CRISPR", fpr); fputs("\n\n", fpr);
            fprintf(fpr, "%s", "FORWARD : 5'(+) -- 3'"); fputs("\n\n", fpr);
            fprintf(fpr, "Pam Length   : %d\n", PAM_LEN);
            fprintf(fpr, "Guide Length : %d\n", GUIDE_LEN);
            fprintf(fpr, "Cut Diff     : %d\n", CUT_DIFF); fputs("\n", fpr);
            fprintf(fpr, "Genome Length : %d\n", seq_len);
            fprintf(fpr, "Total Cut     : %d\n", local_idx); fputs("\n", fpr);
            
            for (int i = 0; i < local_idx; ++i) {
                fprintf(fpr, "Cut-Pos    : %d\n", cut_pos[i]);
                fprintf(fpr, "PAM-Seq    : %.*s\n", PAM_LEN, pam_seq[i]);
                fprintf(fpr, "Target-Seq : %.*s\n", GUIDE_LEN, target_seq[i]);
                fprintf(fpr, "Strand     : %s\n", "forward");
                fputs("\n", fpr);
            }

            fclose(fpr);

            printf("Processor %s, rank %d (master) is finished working\n", processor_name, rank);
        }
        else {
            printf("Can't create or open the file!\n");
            exit(1);
        }

        MPI_Free_mem(cut_pos);
        MPI_Free_mem(pam_seq);
        MPI_Free_mem(target_seq);
    }

    if (rank == SLAVE_NODE)
    {
        MPI_Get_processor_name(processor_name, &name_len);
        
        for (int i = 0; i < seq_len; i++) {

            for (int j = i, idx = 0; j < i + PAM_LEN; j++, idx++) {
                pam_seq_temp[idx] = dna_compliment(*(sequence + j));
            }
            pam_seq_temp[PAM_LEN] = add_Str_Terminator; // add '\0' at the end
                
            // find targets on non-reference (bottom) strand
            // check for out bounds and correct PAM
            if (enough_seq_context(i, seq_len, "REV") && pam_matches(reverse(pam_seq_temp))) {
                // We take the forward of the pam_seq and target_seq to
                // keep the 3` to 5` conventionls

                cut_pos[local_idx] = i + PAM_LEN + CUT_DIFF;

                for (int j = i, pam_sub = 0; j < i + PAM_LEN; j++, pam_sub++) {
                    pam_seq[local_idx][pam_sub] = dna_compliment(*(sequence + j));
                }
                pam_seq[local_idx][PAM_LEN] = add_Str_Terminator;
                    
                for (int j = i + PAM_LEN, tar_sub = 0; j < i + PAM_LEN + GUIDE_LEN; j++, tar_sub++) {
                    target_seq[local_idx][tar_sub] = dna_compliment(*(sequence + j));
                }
                target_seq[local_idx][GUIDE_LEN] = add_Str_Terminator;

                local_idx++;
            }
        }
        
        // write output to file reverse.txt
        FILE *fpr = fopen("reverse.txt", "wb");
        
        if (fpr != NULL) {

            fputs("Programming the genome with CRISPR", fpr); fputs("\n\n", fpr);
            fprintf(fpr, "%s", "REVERSE : 3'(-) -- 5'"); fputs("\n\n", fpr);
            fprintf(fpr, "Pam Length   : %d\n", PAM_LEN);
            fprintf(fpr, "Guide Length : %d\n", GUIDE_LEN);
            fprintf(fpr, "Cut Diff     : %d\n", CUT_DIFF); fputs("\n", fpr);
            fprintf(fpr, "Genome Length : %d\n", seq_len);
            fprintf(fpr, "Total Cut     : %d\n", local_idx); fputs("\n", fpr);
            
            for (int i = 0; i < local_idx; ++i) {
                fprintf(fpr, "Cut-Pos    : %d\n", cut_pos[i]);
                fprintf(fpr, "PAM-Seq    : %.*s\n", PAM_LEN, reverse(pam_seq[i]));
                fprintf(fpr, "Target-Seq : %.*s\n", GUIDE_LEN, reverse(target_seq[i]));
                fprintf(fpr, "Strand     : %s\n", "reverse");
                fputs("\n", fpr);
            }

            fclose(fpr);

            system(SHELLSCRIPT); // shell script to move reverse.txt file from WORKER to MASTER machine.

            printf("Processor %s, rank %d (slave) is finished working\n", processor_name, rank);
        }
        else {

            printf("Can't create or open the file!\n");
            exit(1);
        }

        MPI_Free_mem(cut_pos);
        MPI_Free_mem(pam_seq);
        MPI_Free_mem(target_seq);
    }

    MPI_Free_mem(input_strand);
    MPI_Free_mem(sequence);

    MPI_Barrier(MPI_COMM_WORLD);

    end_time = MPI_Wtime(); // end time

    elapsedTime = end_time - str_time;

    MPI_Reduce(&elapsedTime, &execTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) {
        printf("TOTAL RUNTIME   = %f\n", execTime);
    }

    MPI_Finalize();
}
