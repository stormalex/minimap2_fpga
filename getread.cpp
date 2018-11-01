//#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string>



int main(int argc, char *argv[])
{
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  
  if (argc != 4) {fprintf(stderr,"\nUsage: %s <fasta name> <read name> <read name start line num>\n\n", argv[0]);return -1;}

  char *fasta = argv[1];
  char *readn = argv[2];
  int lineno = atoi(argv[3]);
  if (0 == lineno) lineno = 1;

  fprintf(stderr,"input: fasta [%s], readname [%s], linenum [%d]\n", fasta,readn,lineno);

  fp = fopen(fasta, "r");
  if (fp == NULL) {fprintf(stderr,"Open fasta failed\n"); return -1;}

  int passed = 0;  
  while ((read = getline(&line, &len, fp)) != -1) {
    ++passed;
    if (lineno == passed) {
      printf("%s", line);//readname line
    } else if (lineno < passed) {
      if ('>' == line[0]) break;//next read begin
      printf("%s", line);
    }
  }
  
  free(line);
  return 0;
}
