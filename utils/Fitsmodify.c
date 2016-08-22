// Small program to insert number of rows in fits files to make it close properly

#include<stdio.h>
#include<stdlib.h>
#include"/usr/include/cfitsio/fitsio.h"


int main(int argc, char* argv[]){

// Checking for arguments
if (argc != 3){
  printf("ERROR: invalid arguments\n.");
  printf("USAGE: ./Fitsmodify <FITS file> <number of rows>\n.");
  exit(1);
 }

// Getting the file
fitsfile *fptr;
int status=0,nkeys,ii;
int rows=atoi(argv[2]);
char card[FLEN_CARD];
int data;

//Opening the file
fits_open_file(&fptr, argv[1], READWRITE, &status);
fits_movabs_hdu(fptr, 2, NULL, &status);
fits_update_key(fptr, TINT, "NAXIS2",&rows,NULL,&status);
fits_get_hdrspace(fptr, &nkeys, NULL, &status);
//Printing the header

for (ii = 1; ii <= nkeys; ii++)  {
          fits_read_record(fptr, ii, card, &status); /* read keyword */
          printf("%s\n", card);
     }

//Closing file
fits_close_file(fptr, &status);

if (status)          /* print any error messages */
     fits_report_error(stderr, status);
  return(status);

return(0);


 }

