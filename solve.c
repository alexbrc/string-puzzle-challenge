/*
 * Motivation:
 * http://www.mailund.dk/index.php/2009/03/02/a-string-algorithms-challenge/
 *
 * Expected usage:
 * "Expect the length of x to be from a few hundred thousands to a few millions and k to be between 5 and 20."
 *
 * Target platform:
 * "I will test the programs on a 2Gb RAM dual core i686 Linux box."
 *
 * Fast Karp-Rabin rolling hash function:
 * http://www.eecs.harvard.edu/~ellard/Q-97/PS/sq1997-c8.ps
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct SORTME
{
  int index;
  int count;
};

int compare_sortme(const void *first, const void *second)
{
  return ((struct SORTME *) second)->count - ((struct SORTME *) first)->count;
}

/*
 * The base should be a prime number at least as big as the alphabet
 */
#define BASE 257

/*
 * This function uses the notation from the .ps file mentioned above.
 */
uint64_t get_top_one(int k, int modulus)
{
  uint64_t mask = modulus - 1;
  uint64_t top_one = 1;
  int i;
  for (i=0; i<k; i++)
  {
    top_one = (top_one * BASE) & mask;
  }
  return top_one;
}

/*
 * table: the open-addressed hash table that stores offsets into the data and count array
 * mask: one less than the modulus
 * data: the data read from the file
 * hash_value: the current value of the rolling hash
 * index: the current position in the file
 * count_array: the number of times a pattern has been observed
 * k: the k-mer size
 */
void hash_add(int *table, int mask, uint8_t *data, int hash_value, int index, int *count_array, int k)
{
  int current_position = hash_value;
  int observed_index = table[current_position];
  // move past collisions
  while ( (observed_index != -1) && (memcmp(data + index, data + observed_index, k) != 0) )
  {
    current_position = (current_position + 1) & mask;
    observed_index = table[current_position];
  }
  // initialize the hash entry if necessary
  if (observed_index == -1)
  {
    observed_index = index;
    table[current_position] = observed_index;
  }
  // add a count at the appropriate index
  count_array[observed_index]++;
}

int get_modulus(int file_length)
{
  int m = 1;
  while (m < file_length * 3)
  {
    m <<= 1;
  }
  return m;
}

void do_counts(uint8_t *data, int *count_array, int file_length, int *hash_table, int modulus, int k)
{
  int mask = modulus - 1;
  uint64_t top_one = get_top_one(k, modulus);
  uint64_t hash_value = 0;
  int i;
  // do the initial hash
  for (i=0; i<k; i++)
  {
    hash_value = (hash_value*BASE + data[i]) & mask;
  }
  hash_add(hash_table, mask, data, hash_value, 0, count_array, k);
  // do the rolling hash
  for (i=0; i < file_length - k; i++)
  {
    hash_value = (hash_value*BASE + data[i+k] - data[i]*top_one) & mask;
    hash_add(hash_table, mask, data, hash_value, i+1, count_array, k);
  }
}

int get_nresults(int *count_array, int file_length, int min_count)
{
  // get the number of positions with counts greater than the min count
  int nresults = 0;
  int i;
  for (i=0; i<file_length; i++)
  {
    if (count_array[i] >= min_count)
    {
      nresults += 1;
    }
  }
  return nresults;
}

int main(int argc, char *argv[])
{
  int i;
  // check the number of arguments
  if (argc != 4)
  {
    fprintf(stderr, "Do it like this:\n");
    fprintf(stderr, "%s <filename> <kmer-size> <min-kmer-freq>\n", argv[0]);
    return 1;
  }
  // read the arguments
  char *filename = argv[1];
  int k = atol(argv[2]);
  double min_frequency = atof(argv[3]);
  // check the requested k-mer size
  if (k < 2 || k > 100)
  {
    fprintf(stderr, "expected the k-mer size to be between 2 and 100\n");
    return 1;
  }
  // try to open the file
  FILE *fin = fopen(filename, "rb");
  if (!fin)
  {
    fprintf(stderr, "unable to open file %s\n", filename);
    return 1;
  }
  // get the length of the file
  fseek(fin, 0, SEEK_END);
  int file_length = ftell(fin);
  fseek(fin, 0, SEEK_SET);
  // check the length of the file
  if (file_length > 100000000)
  {
    fprintf(stderr, "the file %s is too long\n", filename);
    return 1;
  }
  // read the file
  char *data = malloc(file_length);
  int nbytes_read = fread(data, file_length, 1, fin);
  fclose(fin);
  // convert the min frequency to a min count
  int max_possible_count = (file_length - k) + 1;
  int min_count = 0;
  if (min_frequency > 0)
  {
    min_count = (int) (min_frequency * max_possible_count + 0.5);
  }
  if (min_count < 1)
  {
    min_count = 1;
  }
  // get an appropriate modulus from the file length
  int modulus = get_modulus(file_length);
  // initialize the count array to all zeros
  int *count_array = calloc(file_length, sizeof(int));
  // initialize the hash table to all invalid indices
  int *hash_table = malloc(modulus * sizeof(int));
  memset(hash_table, 255, modulus * sizeof(int));
  // add the counts
  do_counts(data, count_array, file_length, hash_table, modulus, k);
  // get the number of k-mers that are frequent enough to report
  int nresults = get_nresults(count_array, file_length, min_count);
  // make an results array that will be sorted
  struct SORTME *results_array = malloc(sizeof(struct SORTME) * nresults);
  int current_result_index = 0;
  for (i=0; i<max_possible_count; i++)
  {
    if (count_array[i] >= min_count)
    {
      results_array[current_result_index].count = count_array[i];
      results_array[current_result_index].index = i;
      current_result_index++;
    }
  }
  qsort(results_array, nresults, sizeof(struct SORTME), compare_sortme);
  // print the results
  char *print_buffer = calloc(k+1, 1);
  for (i=0; i<nresults; i++)
  {
    memcpy(print_buffer, data + results_array[i].index, k);
    printf("%s %.12lg\n", print_buffer, results_array[i].count / (double) max_possible_count);
  }
  // free the memory resources
  free(print_buffer);
  free(results_array);
  free(hash_table);
  free(count_array);
  free(data);
  return 0;
}
