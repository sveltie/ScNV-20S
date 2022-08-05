#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "decode.h"
#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_RESET "\x1b[0m"

/*
 *  https://www.ncbi.nlm.nih.gov/nuccore/21557564
 *  LOCUS       NC_004051               2514 bp ss-dna     linear   VRL 13-AUG-2018
 *  DEFINITION  Saccharomyces 20S dna nadnavirus, complete genome.
 *  ACCESSION   NC_004051
 *  VERSION     NC_004051.1
 *  SIZE        2.46 Kb
 */

int main()
{

    FILE *fptr;
    long lSize;
    char *dna_buffer, *comp_buffer, *messenger_buffer;
    char **amino_buffer;
    int i, j, aminoacid_count;
    float buffer_payload;
    size_t base_pairs;
    char *file_name = "ScNV-20S-Genome.txt";

    /* Origin */
    fptr = fopen(file_name, "rb");
    if (!fptr)
        perror(file_name), exit(1);

    fseek(fptr, 0L, SEEK_END);
    lSize = ftell(fptr);
    rewind(fptr);

    /* Allocate memory for entire genome sequence. */
    dna_buffer = calloc(1, lSize + 1);
    if (!dna_buffer)
        fclose(fptr), fputs("memory alloc fails", stderr), exit(1);

    /* Copy file into the buffer */
    if (1 != fread(dna_buffer, lSize, 1, fptr))
        fclose(fptr), free(dna_buffer), fputs("entire read fails", stderr), exit(1);

    for (i = 0; dna_buffer[i] != '\0'; i++)
    {
        /* Remove the number, space, and new line from the genome sequence. */
        while (dna_buffer[i] >= '0' && dna_buffer[i] <= '9' || dna_buffer[i] == ' ' || dna_buffer[i] == '\n')
        {
            for (j = i; dna_buffer[j] != '\0'; j++)
            {
                dna_buffer[j] = toupper(dna_buffer[j + 1]);
            }
            dna_buffer[j] = '\0';
        }
    }

    base_pairs = strlen(dna_buffer);
    aminoacid_count = base_pairs / 3;
    buffer_payload = (float)base_pairs / 1024;
    comp_buffer = (char *)malloc(base_pairs * sizeof(comp_buffer));
    messenger_buffer = (char *)malloc(base_pairs * sizeof(messenger_buffer));
    amino_buffer = NULL;

    /* Read from '5 to '3, coding strand to template strand. */
    invertBuffer(dna_buffer, comp_buffer, base_pairs);
    /* Transcription -> '3 '5 template strand to 5' 3' mRNA. */
    messengerBuffer(comp_buffer, messenger_buffer, base_pairs);
    /* mRNA -> Amino Acid */
    aminoAcidBuffer(messenger_buffer, &amino_buffer, base_pairs);

    printf("Received buffer from\t\t\t%s\n", file_name);
    printf("DNA base pairs(nucleotides)\t\t%ld bp\n", base_pairs);
    printf("Virus genome size\t\t\t%.2f kb\n\n", buffer_payload);
    // printf("ORIGIN (CODING STRAND):\n");
    // printf("%s\n\n", dna_buffer);
    // printf("COMPLEMENTARY DNA (TEMPLATE STRAND):\n");
    // printf("%s\n\n", comp_buffer);
    // printf("mRNA:\n");
    // printf("%s\n\n", messenger_buffer);
    // printf("Protein:\n");
    // printf("%s\n\n", amino_buffer);

    printf("Polypeptide chain or Protein:\n");
    for (size_t i = 0; i < aminoacid_count; i++)
    {
        if (amino_buffer[i] == "M")
            printf(ANSI_COLOR_GREEN "%s" ANSI_COLOR_RESET, amino_buffer[i]);
        else if (amino_buffer[i] == "STOP")
            printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET, amino_buffer[i]);
        else
            printf("%s", amino_buffer[i]);
    }
    printf("\n");

    fclose(fptr);
    free(dna_buffer);
    free(comp_buffer);
    free(messenger_buffer);
    free(amino_buffer);

    return 0;
}

void invertBuffer(char *__buffer, char *__invBuffer, int base_pairs)
{
    /*
     *  Origin -> DNA Complementary
     *
     *  A = Adenine
     *  C = Cytosine
     *  G = Guanine
     *  T = Thymine
     *
     * 5' 3' -> 3' 5'       T -> A
     *                      G <-> C
     */

    int i;
    for (i = 0; i < base_pairs; i++)
    {
        if (__buffer[i] == 'A')
            __invBuffer[i] = 'T';
        else if (__buffer[i] == 'T')
            __invBuffer[i] = 'A';
        else if (__buffer[i] == 'G')
            __invBuffer[i] = 'C';
        else if (__buffer[i] == 'C')
            __invBuffer[i] = 'G';
    }
    __invBuffer[i] = '\0';
}

void messengerBuffer(char *__invBuffer, char *__messengerBuffer, int base_pairs)
{
    /*
     *  Template strand -> Messenger RNA
     *  U = Uracil      T -> A
     *                  G <-> C
     */

    int i;
    for (i = 0; i < base_pairs; i++)
    {
        if (__invBuffer[i] == 'T')
            __messengerBuffer[i] = 'A';
        else if (__invBuffer[i] == 'A')
            __messengerBuffer[i] = 'U';
        else if (__invBuffer[i] == 'G')
            __messengerBuffer[i] = 'C';
        else if (__invBuffer[i] == 'C')
            __messengerBuffer[i] = 'G';
    }
    __messengerBuffer[i] = '\0';
}

int aminoacid_lookup(char __messengerBuffer)
{
    if (__messengerBuffer == 'U')
        return 0;
    else if (__messengerBuffer == 'C')
        return 1;
    else if (__messengerBuffer == 'A')
        return 2;
    else if (__messengerBuffer == 'G')
        return 3;
}

void aminoAcidBuffer(char *__messengerBuffer, char ***__aminoBuffer, int base_pairs)
{
    /*
     *  Ala = A -> GCU, GCC, GCA, GCG
     *  Ile = I -> AUU, AUC, AUA
     *  Arg = R -> CGU, CGC, CGA, CGG; AGA, AGG
     *  Leu = L -> CUU, CUC, CUA, CUG; UUA, UUG
     *  Asn = N -> AAU, AAC
     *  Lys = K -> AAA, AAG
     *  Asp = D -> GAU, GAC
     *  Met = M -> AUG
     *  Phe = F -> UUU, UUC
     *  Cys = C -> UGU, UGC
     *  Pro = P -> CCU, CCC, CCA, CCG
     *  Gln = Q -> CAA, CAG
     *  Ser = S -> UCU, UCC, UCA, UCG; AGU, AGC
     *  Glu = E -> GAA, GAG
     *  Thr = T -> ACU, ACC, ACA, ACG
     *  Trp = W -> UGG
     *  Gly = G -> GGU, GGC, GGA, GGG
     *  Tyr = Y -> UAU, UAC
     *  His = H -> CAU, CAC
     *  Val = V -> GUU, GUC, GUA, GUG
     *  STOP ----> UAA, UGA, UAG
     */

    char *aminoacid_table[4][4][4] =
        {
            {{"F", "F", "L", "L"},
             {"S", "S", "S", "S"},
             {"Y", "Y", "STOP", "STOP"},
             {"C", "C", "STOP", "W"}},
            {{"L", "L", "L", "L"},
             {"P", "P", "P", "P"},
             {"H", "H", "Q", "Q"},
             {"R", "R", "R", "R"}},
            {{"I", "I", "I", "M"},
             {"T", "T", "T", "T"},
             {"N", "N", "K", "K"},
             {"S", "S", "R", "R"}},
            {{"V", "V", "V", "V"},
             {"A", "A", "A", "A"},
             {"D", "D", "E", "E"},
             {"G", "G", "G", "G"}}};

    int i, aa_idx;
    int aminoacid_count = base_pairs / 3;
    *__aminoBuffer = malloc(aminoacid_count * sizeof(**__aminoBuffer));

    for (i = 0, aa_idx = 0; aa_idx < aminoacid_count; i += 3, aa_idx++)
    {
        int x, y, z;
        x = aminoacid_lookup(__messengerBuffer[i]);
        y = aminoacid_lookup(__messengerBuffer[i + 1]);
        z = aminoacid_lookup(__messengerBuffer[i + 2]);
        (*__aminoBuffer)[aa_idx] = aminoacid_table[x][y][z];
    }
}