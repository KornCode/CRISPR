#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
///////////////////
#define arraylen(x) (sizeof (x) / sizeof (*x))
#define newline (printf("\n"))
#define douline (printf("\n\n"))
#define TERMINATOR '\0'

// define constant variables for standard S. Pyogenes Cas9
// for a more flexible program these could be set as function parameters
typedef enum {
    GUIDE_LEN = 20, 
    PAM_LEN = 3, 
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
 Helper function that makes sure the program does not go out of bounds
 trying to find target sequences at the start or end of the sequence.
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
    
    // sequence length, declared as a constant
    const int seq_len = (int) strlen(input_strand);

    double str_time = clock();

    char *sequence = malloc((seq_len + 1) * (sizeof *sequence));
    memcpy(sequence, input_strand, seq_len + 1); // included '\0'

    int max_element = seq_len / PAM_LEN;

    int *cut_pos_fwd = malloc(max_element * (sizeof *cut_pos_fwd));;
    char (*pam_seq_fwd)[PAM_LEN + 1] = malloc((sizeof *pam_seq_fwd) * max_element);
    char (*target_seq_fwd)[GUIDE_LEN + 1] = malloc((sizeof *target_seq_fwd) * max_element);

    int *cut_pos_rev = malloc(max_element * (sizeof *cut_pos_rev));;
    char (*pam_seq_rev)[PAM_LEN + 1] = malloc((sizeof *pam_seq_rev) * max_element);
    char (*target_seq_rev)[GUIDE_LEN + 1] = malloc((sizeof *target_seq_rev) * max_element);
    
    char pam_seq_fwd_temp[PAM_LEN + 1];
    char pam_seq_rev_temp[PAM_LEN + 1];

    int local_idx_fwd = 0;
    int local_idx_rev = 0;

    for (int i = PAM_LEN + GUIDE_LEN; i < seq_len; i++) {
        for (int j = i, idx = 0; j < i + PAM_LEN; j++, idx++) {
            pam_seq_fwd_temp[idx] = *(sequence + j);
        }
        pam_seq_fwd_temp[PAM_LEN] = TERMINATOR; // '\0'

        if (enough_seq_context(i, seq_len, "FWD") && pam_matches(pam_seq_fwd_temp)) {

            cut_pos_fwd[local_idx_fwd] = (i - CUT_DIFF);
                    
            for (int j = i, pam_sub = 0; j < i + PAM_LEN; j++, pam_sub++) {
                pam_seq_fwd[local_idx_fwd][pam_sub] = *(sequence + j);
            }
            pam_seq_fwd[local_idx_fwd][PAM_LEN] = TERMINATOR;
                    
            for (int j = i - GUIDE_LEN, tar_sub = 0; j < i; j++, tar_sub++) {
                target_seq_fwd[local_idx_fwd][tar_sub] = *(sequence + j);
            }
            target_seq_fwd[local_idx_fwd][GUIDE_LEN] = TERMINATOR;

            local_idx_fwd++;
        }
    }

    for (int i = 0; i < seq_len; i++) {
        for (int j = i, idx = 0; j < i + PAM_LEN; j++, idx++) {
            pam_seq_rev_temp[idx] = dna_compliment(*(sequence + j));
        }
        pam_seq_rev_temp[PAM_LEN] = TERMINATOR; // '\0'

        if (enough_seq_context(i, seq_len, "REV") && pam_matches(reverse(pam_seq_rev_temp))) {

            cut_pos_rev[local_idx_rev] = (i + PAM_LEN + CUT_DIFF);

            for (int j = i, pam_sub = 0; j < i + PAM_LEN; j++, pam_sub++) {
                pam_seq_rev[local_idx_rev][pam_sub] = dna_compliment(*(sequence + j));
            }
            pam_seq_rev[local_idx_rev][PAM_LEN] = TERMINATOR;
                    
            for (int j = i + PAM_LEN, tar_sub = 0; j < i + PAM_LEN + GUIDE_LEN; j++, tar_sub++) {
                target_seq_rev[local_idx_rev][tar_sub] = dna_compliment(*(sequence + j));
            }
            target_seq_rev[local_idx_rev][GUIDE_LEN] = TERMINATOR;

            local_idx_rev++;
        }
    }

    FILE *fpf = fopen("forward_seq.txt", "wb");

    if (fpf != NULL) {
        fputs("Programming the genome with CRISPR", fpf); fputs("\n\n", fpf);
        fprintf(fpf, "%s", "FORWARD : 5'(+) -- 3'"); fputs("\n\n", fpf);
        fprintf(fpf, "Pam Length   : %d\n", PAM_LEN);
        fprintf(fpf, "Guide Length : %d\n", GUIDE_LEN);
        fprintf(fpf, "Cut Diff     : %d\n", CUT_DIFF); fputs("\n", fpf);
        fprintf(fpf, "Genome Length : %d\n", seq_len);
        fprintf(fpf, "Total Cut     : %d\n", local_idx_fwd); fputs("\n", fpf);

        for (int i = 0; i < local_idx_fwd; ++i) {
            fprintf(fpf, "Cut-Pos    : %d\n", cut_pos_fwd[i]);
            fprintf(fpf, "PAM-Seq    : %.*s\n", PAM_LEN, pam_seq_fwd[i]);
            fprintf(fpf, "Target-Seq : %.*s\n", GUIDE_LEN, target_seq_fwd[i]);
            fprintf(fpf, "Strand     : %s\n", "forward");
            fputs("\n", fpf);
        }

        fclose(fpf);
    }
    else {
        printf("Can't create or open the file!\n");
        exit(1);
    }

    FILE *fpr = fopen("reverse_seq.txt", "wb");

    if (fpr != NULL)
    {
        fputs("Programming the genome with CRISPR", fpr); fputs("\n\n", fpr);
        fprintf(fpr, "%s", "REVERSE : 3'(-) -- 5'"); fputs("\n\n", fpr);
        fprintf(fpr, "Pam Length   : %d\n", PAM_LEN);
        fprintf(fpr, "Guide Length : %d\n", GUIDE_LEN);
        fprintf(fpr, "Cut Diff     : %d\n", CUT_DIFF); fputs("\n", fpr);
        fprintf(fpr, "Genome Length : %d\n", seq_len);
        fprintf(fpr, "Total Cut     : %d\n", local_idx_rev); fputs("\n", fpr);

        for (int i = 0; i < local_idx_rev; ++i)
        {
            fprintf(fpr, "Cut-Pos    : %d\n", cut_pos_rev[i]);
            fprintf(fpr, "PAM-Seq    : %.*s\n", PAM_LEN, reverse(pam_seq_rev[i]));
            fprintf(fpr, "Target-Seq : %.*s\n", GUIDE_LEN, reverse(target_seq_rev[i]));
            fprintf(fpr, "Strand     : %s\n", "reverse");
            fputs("\n", fpr);
        }

        fclose(fpr);
    }
    else {
        printf("Can't create or open the file!\n");
        exit(1);
    }

    free(input_strand); 
    free(sequence);

    free(cut_pos_fwd);
    free(pam_seq_fwd);
    free(target_seq_fwd);
    
    free(cut_pos_rev);
    free(pam_seq_rev);
    free(target_seq_rev);

    double end_time = clock();

    printf("%f\n", (end_time - str_time) / CLOCKS_PER_SEC);
}







