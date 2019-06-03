#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "numrec/include/nr.h"
#include <cfitsio/fitsio.h>
#include <cfitsio/fitsio2.h>
#include "/home/arasmus/ki-rhedwig/arasmus/MIT/new/rv.h"

typedef enum {_exclude,_bias,_image} filtertype;

typedef struct bias_array {
  int b_nelem;
  float *b_lvl;
} bias_array;
  
typedef struct img_basepars {
  float bias; // bias will be either the avg of specified overscan pixels or avg of pixels in upper right region in the case of bias_model
  float gain;
  bias_array model_bias_array[2];
  float model_bias;
} img_basepars;

typedef struct extdata_map {
  int data_swap_hdu[64]; // data_swap will be populated later: pairs of indices [2*ix+0],[2*ix+1]
  char **swap_map;      // pairs of EXTNAME keywords for matching: indices [2*ix+0],[2*ix+1]
  int nswap;
  int data_ext_map[32];
} extdata_map;

typedef struct img_stack {
  fitsfile **ffp;
  char     **ffname;
  fitsfile *bias_ffp;
  char     *bias_ffname;
  int      *hdu;
  img_basepars *ibp;
  char *output_file;
  fitsfile *output_ffp;
  FILE *diagnostic_output;
  FILE *histogram_output;
  FILE *img_basepar_output, *img_basepar_input;
  FILE *xray_output_file;
  char *xray_output_file_prefix;
  float *bias_by_hdu, *scale_by_hdu;
  char *incl;
  int bias_subtract;
  int bias_1d_plus_1d_model;
  int bias_min_ser, bias_min_par;
  int ser_oscan_min,ser_oscan_max,par_oscan_min,par_oscan_max;
  int bias_model_smooth_npix;
  int data_min_par, data_max_par, data_min_ser, data_max_ser;
  int normalize;
  int integer;
  int Xray;
  int segment_by_segment,multiplex_mosaic;
  int verbose;
  int mult,diff,sum,average,median;
  int n_ff_chunk;
  int n_ff;
  int fits_status;
  extdata_map edm;
} img_stack;

void usage (void);
void printerror( int status);
int  float_compare(const void *a, const void *b);
void compute_image_base_pars(img_stack *isk,int im_ix);
void get_ibp_list(img_stack *isk);
void img_isr(img_stack *isk,int ff_ix,float *im,long *naxes);
void img_isr_chunk(img_stack *isk,int ff_ix,float *im,long *naxes,
		   long *fpixel_read,long npixel_read);

// fitsfile *ff,int minser,img_basepars *ibp);

char *ccd250_088_swaplist[]={"Segment00","Segment07",
			     "Segment01","Segment06",
			     "Segment02","Segment05",
			     "Segment03","Segment04"};

int main(int argc,char **argv) {

  static img_stack is;

  is.verbose=0;
  is.edm.nswap=0;
  is.incl=NULL;
  is.diagnostic_output=NULL;
  is.img_basepar_output=NULL;
  is.img_basepar_input=NULL;
  is.histogram_output=NULL;
  is.bias_by_hdu=NULL;
  is.scale_by_hdu=NULL;
  // for now expect a list of input files on the command line.
  is.n_ff_chunk=4;
  is.n_ff=0;  
  is.integer=0;  
  is.normalize=0;  
  is.Xray=0;
  is.xray_output_file_prefix=NULL;
  is.segment_by_segment=1;
  is.multiplex_mosaic=0;  

  is.bias_subtract=0;
  is.bias_1d_plus_1d_model=0;

  is.bias_min_ser=5;
  is.bias_min_par=5;

  is.data_min_par=1;
  is.data_max_par=50000;
  is.data_min_ser=1;
  is.data_max_ser=50000;

  is.mult = is.diff = is.sum = is.average = is.median = 0;
  is.output_file=NULL;
  is.fits_status=0;
  is.ffp=(fitsfile**)malloc(is.n_ff_chunk*sizeof(fitsfile*));
  is.ffname=(char**)malloc(is.n_ff_chunk*sizeof(char*));

  is.bias_ffname=NULL;
  is.bias_ffp=NULL;

  is.hdu=(int*)malloc(is.n_ff_chunk*sizeof(int));
  is.ibp=(img_basepars*)malloc(is.n_ff_chunk*sizeof(img_basepars));
  
  extern int optind,opterr,optopt;
  extern char *optarg;
  optind=1;
  int c;

  static struct option long_options[] = {
    {"verbose",           no_argument, NULL, 'v'},
    {"mult",              no_argument, NULL, 'p'},
    {"diff",              no_argument, NULL, 'D'},
    {"sum",               no_argument, NULL, 's'},
    {"average",           no_argument, NULL, 'a'},
    {"median",            no_argument, NULL, 'm'},
    {"integer",           no_argument, NULL, 'i'},
    {"normalize",         no_argument, NULL, 'z'},
    {"Xray",              no_argument, NULL, 'X'},
    {"multiplex",         no_argument, NULL, 'M'},
    {"bias_subtract",     no_argument, NULL, 'b'},
    {"bias_model",        no_argument, NULL, 'd'},
    {"bias_min_parallel", required_argument, NULL, 0x09},
    {"bias_min_serial",   required_argument, NULL, 0x0A},
    {"data_min_parallel", required_argument, NULL, 0x0B},
    {"data_max_parallel", required_argument, NULL, 0x0C},
    {"data_min_serial",   required_argument, NULL, 0x0D},
    {"data_max_serial",   required_argument, NULL, 0x0E},
    {"diagnostic_output", required_argument, NULL, 0x0F},
    {"img_basepar_output",required_argument, NULL, 0x10},
    {"img_basepar_input", required_argument, NULL, 0x11},
    {"88",                no_argument,       NULL, 0x12}, // correct mapping of CCD250-088
    {"Xray_output_file_prefix", required_argument, NULL, 0x13},
    {"bias_file",   required_argument,       NULL, 0x14},
    {"output_file", required_argument,       NULL, 'o'},
    {"histogram_output", required_argument,  NULL, 'h'},
    {0,0,0,0}
  };

  while ((c=getopt_long(argc,argv,"dvzXMbsamibo:",long_options,NULL))!=-1) {
    switch(c) {
    case 0:
      // switch assigned. keep going.
      break;
    case 'v':
      is.verbose=1;
      break;
    case 'z':
      is.normalize=1;
      break;
    case 'X':
      is.Xray=1;
      break;
    case 'M':
      // toggle segment-by-segment & multiplex_mosaic
      is.multiplex_mosaic=1;
      is.segment_by_segment=0;
      break;
    case 'd':
      is.bias_1d_plus_1d_model=1;
      break;
    case 'b':
      is.bias_subtract=1;
      break;
    case 'p':
      is.mult=1;
      break;
    case 'D':
      is.diff=1;
      break;
    case 's':
      is.sum=1;
      break;
    case 'a':
      is.average=1;
      break;
    case 'm':
      is.median=1;
      break;
    case 'i':
      is.integer=1;
      break;
    case 0x09:
      is.bias_min_par=atoi(optarg);
      break;
    case 0x0A:
      is.bias_min_ser=atoi(optarg);
      break;
    case 0x0B:
      is.data_min_par=atoi(optarg);
      break;
    case 0x0C:
      is.data_max_par=atoi(optarg);
      break;
    case 0x0D:
      is.data_min_ser=atoi(optarg);
      break;
    case 0x0E:
      is.data_max_ser=atoi(optarg);
      break;
    case 0x0F:
      if ((is.diagnostic_output=fopen(optarg,"w"))==NULL) {
	fprintf(stderr,"can't open diagnostic file %s.\n",optarg);
	exit(1);
      }
      break;
    case 0x10:
      if ((is.img_basepar_output=fopen(optarg,"w"))==NULL) {
	fprintf(stderr,"can't open img_basepar_output file %s.\n",optarg);
	exit(1);
      }
      break;
    case 0x11:
      if ((is.img_basepar_input=fopen(optarg,"r"))==NULL) {
	fprintf(stderr,"can't open img_basepar_input file %s.\n",optarg);
	exit(1);
      }
      break;
    case 0x12:
      // use swap mapping of char *ccd250_088_swaplist[].
      is.edm.swap_map=ccd250_088_swaplist;
      if (sizeof(is.edm.swap_map) % 2 != 0) {
	fprintf(stderr,"odd number of entries in swaplist:\n");
	fprintf(stderr,"%ld\n",sizeof(is.edm.swap_map));
	exit(1);
      }
      is.edm.nswap = sizeof(is.edm.swap_map);
      break;
    case 0x13:
      is.xray_output_file_prefix=optarg;
      break;
    case 0x14:
      is.bias_ffname=optarg;
      break;
    case 'o':
      is.output_file=optarg;
      break;
    case 'h':
      if ((is.histogram_output=fopen(optarg,"w"))==NULL) {
	fprintf(stderr,"can't open histogram output file %s.\n",optarg);
	exit(1);
      }
      break;
    case 1:
      fprintf(stderr,"regular argument: %s\n",optarg);
      break;
    case ':':
    case '?':
      // getopt_long already printed an error message..
      fprintf(stderr,"getopt returns a problem: %02x (%c) (optarg=%s)\n",c,c,optarg);
      usage();
      break;
    default:
      abort();
    }
  }

  if (is.Xray) {
    // override some of the modes.
    is.median=1;
    if (is.average)
      fprintf(stderr,"Xray mode overrides --average. resetting..\n");
    if (is.sum)
      fprintf(stderr,"Xray mode overrides --sum. resetting..\n");
    if (!is.median)
      fprintf(stderr,"Xray mode implies --median. setting..\n");
    if (is.multiplex_mosaic)
      fprintf(stderr,"Xray mode overrides --multiplex_mosaic. resetting..\n");
    if (!is.bias_subtract)
      fprintf(stderr,"Xray mode implies --bias_subtract. setting..\n");
    if (is.normalize)
      fprintf(stderr,"Xray mode overrides --normalize. resetting..\n");
    is.average = is.sum = 0;
    is.multiplex_mosaic = 0;
    is.bias_subtract = 1;
    is.normalize = 0;
  }

  // some higher level checks on run parameters
  if ((is.sum + is.average + is.median + is.diff + is.mult)>1) {
    fprintf(stderr,"specify only ONE of {sum, average, median, diff, mult}\n");
    exit(1);
  }

  if ((is.sum + is.average + is.median + is.diff + is.mult
       + is.img_basepar_output)<=0) {
    fprintf(stderr,"nothing to do! need to specify one of {sum, average, median, diff}\n");
    exit(1);
  }

  if ((is.bias_subtract + is.bias_1d_plus_1d_model > 0) &&
      (is.bias_subtract * is.bias_1d_plus_1d_model != 0)) {
    fprintf(stderr,"optionally specify either (not both) of {bias_subtract, bias_model}\n");
    exit(1);
  }

  if ((is.output_file==NULL) && !is.img_basepar_output) {
    fprintf(stderr,"no output file set! specify with --output_file=outfile\n");
    //    exit(1);
  }

  if (is.verbose)
    puts("verbose flag is set");
  if (optind<argc) {
    printf ("non-option ARGV-elements: ");
    while (optind<argc) {
      printf ("%s ",argv[optind]);
      // append filename list:
      if (is.n_ff >= is.n_ff_chunk) {
	// grow lists as necessary
	is.n_ff_chunk *= 2;
	is.ffname=(char**)realloc(is.ffname,is.n_ff_chunk*sizeof(char*));
	is.ffp=(fitsfile**)realloc(is.ffp,is.n_ff_chunk*sizeof(fitsfile*));
	is.hdu=(int*)realloc(is.hdu,is.n_ff_chunk*sizeof(int));
	is.ibp=(img_basepars*)realloc(is.ibp,is.n_ff_chunk*sizeof(img_basepars));
      }
      is.ffname[is.n_ff]=argv[optind];
      is.n_ff++;
      optind++;
    }
    putchar('\n');
  }

  // read in img_basepar_input values, if available. then close the file.
  if (is.img_basepar_input) get_ibp_list(&is);
  

  // now go through list to make sure they were registered.
  int i;
  fprintf(stderr,"my list of %d files:\n",is.n_ff);
  
  int goodfiles=0;
  int badfiles=0;
  // go through the source files
  for (i=0;i<is.n_ff;i++) {
    fprintf(stderr,"%s ",is.ffname[i]);
    // open and populate the fitsfile object for each filename
    if (fits_open_file(&is.ffp[i],is.ffname[i],READONLY,&is.fits_status)) 
      badfiles++;
    else 
      goodfiles++;
    is.ibp[i].bias=0.0;
    is.ibp[i].gain=1.0;
  }
  // and check the specified bias file, if any
  if (is.bias_ffname) {
    fprintf(stderr,"\nbias input file specified: %s\n",is.bias_ffname);
    if (fits_open_file(&is.bias_ffp,is.bias_ffname,READONLY,&is.fits_status)) 
      badfiles++;
    else
      goodfiles++;
  }

  if (badfiles) {
    fprintf(stderr,"appear to have %d good files and %d bad files.\n",
	    goodfiles,badfiles);
    printerror(is.fits_status);
  } else {
    fprintf(stderr,"opened %d files.\n",goodfiles);
  }

  // if operation is is.diff=1 make sure that goodfiles is equal to 2.
  // else exit.

  if (is.diff || is.mult) {
    if (goodfiles !=2 ) {
      if (is.diff)
	fprintf(stderr,"\nasking for difference but cmdline includes %d files. "
		"what to do!?\n\n",goodfiles);
      if (is.mult)
	fprintf(stderr,"\nasking for product but cmdline includes %d files. "
		"what to do!?\n\n",goodfiles);
      usage();
    } else {
      fprintf(stderr,"got %d files to difference\n",goodfiles);
    }
  }
  if (is.fits_status)
    printerror(is.fits_status);
  // so far so good.
  if (1) { // report the structure of the fits file and exit.
    int nhdu=0;
    if (fits_get_num_hdus(is.ffp[0],&nhdu,&is.fits_status))
      printerror(is.fits_status);
    int hdu_ix;
    int hdu_type;
    int morekeys=0;
    for (hdu_ix=1;hdu_ix<=nhdu;hdu_ix++) {
      if (fits_movabs_hdu(is.ffp[0],hdu_ix,&hdu_type,&is.fits_status))
	printerror(is.fits_status);
    }
    //    exit(0);
  }
  
  // now open the output file
  if (is.output_file) {
    fprintf(stderr,"outputfile = %s\n",is.output_file);
    if (fits_create_file(&is.output_ffp,is.output_file,&is.fits_status)) 
      printerror(is.fits_status);
    
    // as an example, copy entire file from first input to the output.
    int nhdu=0;
    if (fits_get_num_hdus(is.ffp[0],&nhdu,&is.fits_status))
      printerror(is.fits_status);
    int hdu_ix;
    int hdu_type;
    int morekeys=0;

    // go through 1st input file once to determine the mapping indices.
    for (hdu_ix=1;hdu_ix<=nhdu;hdu_ix++) {

      if (fits_movabs_hdu(is.ffp[0],hdu_ix,&hdu_type,&is.fits_status))
	printerror(is.fits_status);

      is.edm.data_ext_map[hdu_ix]=hdu_ix;

      if (hdu_type==IMAGE_HDU) {
	char extname[2048];
	if (fits_read_key(is.ffp[0],TSTRING,"EXTNAME",
			  extname,NULL,&is.fits_status)) {
	  if (is.fits_status != KEY_NO_EXIST) 
	    printerror(is.fits_status);
	  else
	    is.fits_status=0;
	} else {	  
	  // go through the list of swaps, if any, and populate.
	  int swap_ix=is.edm.nswap;
	  while (swap_ix--) {
	    if (strcmp(is.edm.swap_map[swap_ix],extname)==0) {
	      // latch this hdu
	      is.edm.data_swap_hdu[swap_ix]=hdu_ix;
	    }
	  }
	}
      }
    }

    // first pass is finished.
    {
      // now look for swap_map pairs to properly populate is.edm.data_ext_map[]
      int swap_ix;
      for (swap_ix=0;swap_ix<is.edm.nswap;swap_ix+=2) {
	int criss=is.edm.data_swap_hdu[swap_ix];
	int cross=is.edm.data_swap_hdu[swap_ix+1];
	is.edm.data_ext_map[criss]=cross;
	is.edm.data_ext_map[cross]=criss;
      }
      fprintf(stderr,"HDU mapping (for image data only):\n");
      for (hdu_ix=1;hdu_ix<=nhdu;hdu_ix++) {
	fprintf(stderr," %d",is.edm.data_ext_map[hdu_ix]);
      }
      fprintf(stderr,"\n");
    }

    int *shipped_out_hdu=(int*)malloc((nhdu+1)*sizeof(int)); // used only when combining input & output img_basepars
    
    for (hdu_ix=1;hdu_ix<=nhdu;hdu_ix++) {
      shipped_out_hdu[hdu_ix]=0;
      int ff_ix;
      for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
	if (fits_movabs_hdu(is.ffp[ff_ix],hdu_ix,&hdu_type,&is.fits_status))
	  printerror(is.fits_status);
	is.hdu[ff_ix]=hdu_ix;
      }

      // get HDU type from is.ffp[0]
      fits_get_hdu_type(is.ffp[0],&hdu_type,&is.fits_status);
      // previously is.output_ffp had been copied from is.ffp[0]
      //      fits_get_hdu_type(is.output_ffp,&hdu_type,&is.fits_status);

      long totpix=1L;
      int naxis;
      int bitpix;
      int maxdim=9;
      long naxes[maxdim];

      if (hdu_type == IMAGE_HDU) {
	// compute the # of pixels
	int ii;
	for (int ii=0;ii<maxdim;ii++) naxes[ii]=1L;
	fits_get_img_param(is.ffp[0],maxdim,&bitpix,&naxis,naxes,&is.fits_status);
	for (ii=0;ii<maxdim;ii++)
	  totpix *= naxes[ii];
      }
      
      if ((hdu_type != IMAGE_HDU) || (naxis==0) || (totpix==0)) {
	// copy tables & null images
	if (fits_copy_hdu(is.ffp[0],
			  is.output_ffp,morekeys,&is.fits_status))
	  printerror(is.fits_status);
      } else {
	// hdu_type==IMAGE_HDU and contains data. need to explicitly
	// create a new image supporting compression
	if (fits_create_img(is.output_ffp,bitpix,naxis,naxes,&is.fits_status))
	  printerror(is.fits_status);

	// first test whether this is a compressed image. If not, 
	if (fits_is_compressed_image(is.output_ffp,&is.fits_status)) {
	  // image is not compressed
	  int tstatus=0;
	  char card[81];
	  fits_read_card(is.ffp[0],"EXTNAME",card,&tstatus);
	  if (tstatus) {
	    strcpy(card,
	   "EXTNAME = 'COMPRESSED_IMAGE'   / name of this binary table extension");
	    if (fits_write_record(is.output_ffp,card,&is.fits_status))
	      printerror(is.fits_status);
	  }
	  if (0) {
	    int comptype;
	    fprintf(stderr,"this image is compressed. type = %d\n",
		    fits_get_compression_type(is.output_ffp,&comptype,&is.fits_status));
	  }
	} else {
	  // image is not compressed
	}
	
	// copy non-structural keywords into the output
	int nkeys,ii;
	fits_get_hdrspace(is.ffp[0],&nkeys,NULL,&is.fits_status);
	for (ii=1;ii<=nkeys;ii++) {
	  char card[81];
	  fits_read_record(is.ffp[0],ii,card,&is.fits_status);
	  if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	    fits_write_record(is.output_ffp,card,&is.fits_status);
	}
	// so far we're probably close to imcopy.c line 182..
	// reverse the values and write back into the file
	int maxdim=5;
	int bitpix,naxis;
	long naxes[maxdim];
	fits_get_img_param(is.output_ffp,maxdim,
			   &bitpix,&naxis,naxes,&is.fits_status);
	if (naxis!=2) continue;
	if (naxes[0]==0 || naxes[1]==0) continue;
	// check for EXTNAME
	char extname[2048];
	int  chan_num;
	if (fits_read_key(is.output_ffp,TSTRING,"EXTNAME",
			  extname,NULL,&is.fits_status)) {
	  if (is.fits_status==KEY_NO_EXIST)
	    is.fits_status=0;
	  else 
	    printerror(is.fits_status);
	}

	if (strcmp(extname,"MOSAIC")==0) continue;

	if (fits_read_key(is.output_ffp,TINT,"CHANNEL",
			  &chan_num,NULL,&is.fits_status)) {
	  fprintf(stderr,"channel maybe not defined? chan_num = %d\n",chan_num);
	  // try keyword IMAGEID - used by ITL
	  is.fits_status=0;
	  if (fits_read_key(is.output_ffp,TINT,"IMAGEID",
			    &chan_num,NULL,&is.fits_status)) {
	    fprintf(stderr,"imageid doesn't work either: chan_num = %d\n",chan_num);
	    is.fits_status=0;
	    if (fits_read_key(is.output_ffp,TINT,"AMPNO",
			      &chan_num,NULL,&is.fits_status)) {
	      fprintf(stderr,"ampno doesn't work either: chan_num = %d\n",chan_num);
	      exit(0);
	    }
	  } 
	  fprintf(stderr,"will use chan_num = %d\n",chan_num);
	}

	//	fits_set_bscale(is.output_ffp,2,0,&is.fits_status);
	// read, flip, rewrite
	float *im,*output_im;
	long fpixel[]={1,1};
	im=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
	output_im=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));

	if ((im==NULL) || (output_im==NULL)) {
	  fprintf(stderr,"can't allocate!\n");
	  usage();
	}

	long m; // pixel index for use below

	// input file trumps any img_base_par computation
	// but if is.bias_ffname is specified, the pixel-by-pixel
	// values will trump is.img_basepar_input

	for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
	  if (is.bias_subtract || is.bias_1d_plus_1d_model || is.img_basepar_output) 
	    compute_image_base_pars(&is,ff_ix);
	  if (is.img_basepar_input) {
	    is.ibp[ff_ix].bias = is.bias_by_hdu[is.hdu[ff_ix]];
	    is.ibp[ff_ix].gain = is.scale_by_hdu[is.hdu[ff_ix]];
	    if (is.img_basepar_output && !shipped_out_hdu[is.hdu[ff_ix]]) {
	      // output the effective biases & scales used, which may be useful for tweaking
	      fprintf(is.img_basepar_output,"filename %s[%s] HDU %d bias %g scale %g\n",
		      "BASED_ON_INPUT_IBP_FILE",extname,
		      is.hdu[ff_ix],is.ibp[ff_ix].bias,is.ibp[ff_ix].gain);
	      fflush(is.img_basepar_output);
	      shipped_out_hdu[is.hdu[ff_ix]]=1;
	    }
	  }
	}
	
	if (is.Xray) {
	  char output_xray_file[4096];
	  if (is.xray_output_file_prefix) {
	    sprintf(output_xray_file,"%s_%s_CH%d.evlist",
		    is.xray_output_file_prefix,extname,chan_num);
	  } else {
	    sprintf(output_xray_file,"%s_CH%d.evlist",extname,chan_num);
	  }
	  is.xray_output_file=fopen(output_xray_file,"w");
	  if (is.xray_output_file==NULL) {
	    fprintf(stderr,"can't open xray output file.\n");
	    usage();
	  }
	}

	if (is.average || is.sum) {
	  // zero output image:
	  m=naxes[0]*naxes[1];
	  while (m--) output_im[m]=0.0;
	  // cycle through input files here to produce resulting image content:
	  for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
	    {
	      int hdu_save;
	      fits_get_hdu_num(is.ffp[ff_ix],&hdu_save);
	      if (fits_movabs_hdu(is.ffp[ff_ix],is.edm.data_ext_map[hdu_save],
				  NULL,&is.fits_status) || 
		  fits_read_pix(is.ffp[ff_ix],TFLOAT,fpixel,
				(long)naxes[0]*naxes[1],
				NULL,im,NULL,&is.fits_status) || 
		  fits_movabs_hdu(is.ffp[ff_ix],hdu_save,
				  NULL,&is.fits_status)) 
		printerror(is.fits_status);
	      // probably not so efficient if working with a stack of input
	      // images, but pixel-by-pixel bias subtraction can work
	      // here if the is.bias_ffname was specified.
	      if (is.bias_ffname) {
		float *bias_im=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
		if (bias_im==NULL) {
		  fprintf(stderr,"can't allocate for bias_im.\n");
		  exit(1);
		}
		if (fits_movabs_hdu(is.bias_ffp,is.edm.data_ext_map[hdu_save],
				    NULL,&is.fits_status) || 
		    fits_read_pix(is.bias_ffp,TFLOAT,fpixel,
				  (long)naxes[0]*naxes[1],
				  NULL,bias_im,NULL,&is.fits_status) || 
		    fits_movabs_hdu(is.bias_ffp,hdu_save,
				    NULL,&is.fits_status)) 
		  printerror(is.fits_status);
		// and now subtract bias_im from im, pixel by pixel.
		long pix_ix=naxes[0]*naxes[1];
		while (pix_ix--)
		  im[pix_ix] -= bias_im[pix_ix];
		free(bias_im);
	      }
	    }

	    // subtract bias & scale image just read
	    if ((is.bias_1d_plus_1d_model != 0) ||
		(is.ibp[ff_ix].bias != 0) ||
		(is.ibp[ff_ix].gain != 1.0)) 

	      img_isr(&is,ff_ix,im,naxes);

	    // increment the output image
	    m=naxes[0]*naxes[1];
	    while (m--) 
	      output_im[m] += im[m];
	  }
	  // if average image is requested, divide by the number of input images
	  if (is.average || is.integer) {
	    m=naxes[0]*naxes[1];
	    while (m--) {
	      if (is.average)
		output_im[m] /= is.n_ff;
	      if (is.integer)
		output_im[m] = floor(output_im[m]+0.5);
	    }
	  }
	}

	if (is.diff || is.mult) {
	  // zero output image:
	  m=naxes[0]*naxes[1];
	  while (m--) output_im[m]=0.0;
	  // cycle through input files here to produce resulting image content:
	  for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
	    {
	      int hdu_save;
	      fits_get_hdu_num(is.ffp[ff_ix],&hdu_save);
	      if (fits_movabs_hdu(is.ffp[ff_ix],is.edm.data_ext_map[hdu_save],
				  NULL,&is.fits_status) || 
		  fits_read_pix(is.ffp[ff_ix],TFLOAT,fpixel,
				(long)naxes[0]*naxes[1],
				NULL,im,NULL,&is.fits_status) || 
		  fits_movabs_hdu(is.ffp[ff_ix],hdu_save,
				  NULL,&is.fits_status)) 
		printerror(is.fits_status);
	      // probably not so efficient if working with a stack of input
	      // images, but pixel-by-pixel bias subtraction can work
	      // here if the is.bias_ffname was specified.
	      if (is.bias_ffname) {
		float *bias_im=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
		if (bias_im==NULL) {
		  fprintf(stderr,"can't allocate for bias_im.\n");
		  exit(1);
		}
		if (fits_movabs_hdu(is.bias_ffp,is.edm.data_ext_map[hdu_save],
				    NULL,&is.fits_status) || 
		    fits_read_pix(is.bias_ffp,TFLOAT,fpixel,
				  (long)naxes[0]*naxes[1],
				  NULL,bias_im,NULL,&is.fits_status) || 
		    fits_movabs_hdu(is.bias_ffp,hdu_save,
				    NULL,&is.fits_status)) 
		  printerror(is.fits_status);
		// and now subtract bias_im from im, pixel by pixel.
		long pix_ix=naxes[0]*naxes[1];
		while (pix_ix--)
		  im[pix_ix] -= bias_im[pix_ix];
		free(bias_im);
	      }
	    }

	    // subtract bias & scale image just read
	    if ((is.bias_1d_plus_1d_model != 0) ||
		(is.ibp[ff_ix].bias != 0) ||
		(is.ibp[ff_ix].gain != 1.0)) 

	      img_isr(&is,ff_ix,im,naxes);

	    // increment the output image
	    m=naxes[0]*naxes[1];
	    if (is.diff) { // difference
	      if (ff_ix%2==0) 
		while (m--) output_im[m] += im[m];
	      else
		while (m--) output_im[m] -= im[m];
	    }
	    if (is.mult) { // product
	      if (ff_ix%2==0) 
		while (m--) output_im[m] = im[m];
	      else
		while (m--) output_im[m] *= im[m];
	    }
	  }
	  // if integer output format is requested, round result
	  if (is.integer) {
	    m=naxes[0]*naxes[1];
	    while (m--) {
	      if (is.integer)
		output_im[m] = floor(output_im[m]+0.5);
	    }
	  }
	}

	if (is.median) {
	  // not sure how buffered the CFITSIO system is in regards to pixel data.
	  // 1 pixel at a time
	  long m_max=naxes[0]*naxes[1];
	  float *median_sample=(float*)malloc(is.n_ff*sizeof(float));


	  if (median_sample==NULL) {
	    fprintf(stderr,"can't allocate.\n");
	    usage();
	  }


	  // compute the median by reading in a chunk of image at a time
	  // from each input file
	  long mychunk=m_max;
	  mychunk=8192;
	  // declare buffers for each image
	  float **imbuf=NULL;
	  imbuf=(float**)malloc(is.n_ff*sizeof(float*));
	  if (imbuf==NULL) {
	    fprintf(stderr,"can't allocate.\n");
	    usage();
	  }
	  int buf_ix;
	  for (buf_ix=0;buf_ix<is.n_ff;buf_ix++) {
	    if ((imbuf[buf_ix]=(float*)malloc(mychunk*sizeof(float)))==NULL) {
	      fprintf(stderr,"can't allocate.\n");
	      usage();
	    }
	  }
	  long thism;
	  for (thism=0;thism<m_max;thism+=mychunk) {
	    // read in segment from each file
	    long fpixel_read[]={thism%naxes[0]+1,thism/naxes[0]+1};
	    long npixel_read=((naxes[0]*naxes[1]-thism>mychunk)?
			      mychunk:
			      (naxes[0]*naxes[1]-thism));
	    // read from files
	    for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {

	      {
		int hdu_save;
		fits_get_hdu_num(is.ffp[ff_ix],&hdu_save);
		if (fits_movabs_hdu(is.ffp[ff_ix],is.edm.data_ext_map[hdu_save],
				    NULL,&is.fits_status) ||
		    fits_read_pix(is.ffp[ff_ix],TFLOAT,
				  fpixel_read,npixel_read,
				  NULL,imbuf[ff_ix],NULL,
				  &is.fits_status) ||
		    fits_movabs_hdu(is.ffp[ff_ix],hdu_save,
				    NULL,&is.fits_status))
		  printerror(is.fits_status);
	      }

	      // subtract bias & scale image just read
	      if ((is.bias_1d_plus_1d_model != 0) ||
		  (is.ibp[ff_ix].bias != 0) ||
		  (is.ibp[ff_ix].gain != 1.0)) {
		img_isr_chunk(&is,ff_ix,imbuf[ff_ix],naxes,fpixel_read,npixel_read);
	      }
	    }

	    for (m=thism;m<thism+npixel_read;m++) {
	      for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
		median_sample[ff_ix]=imbuf[ff_ix][m-thism];
	      }
	      // sort the list and choose the central value
	      qsort(median_sample,is.n_ff,sizeof(float),float_compare);
	      if (is.n_ff%2==1) {
		output_im[m]=
		  0.5*(median_sample[is.n_ff/2]+median_sample[is.n_ff/2+1]);
	      } else {
		output_im[m]=median_sample[is.n_ff/2];
	      }
	    }
	  }
	  // deallocate.
	  for (buf_ix=0;buf_ix<is.n_ff;buf_ix++) 
	    free(imbuf[buf_ix]);
	  free(imbuf);
	  free(median_sample);
	}

	// prepare & write the result..

	// resize the output image. this may be necessary to avoid overflow..
	if (0) {
	  float a;
	  fits_get_version(&a);
	  fprintf(stderr,"fits version is %f\n",a);
	  exit(0);
	}

	// this was throwing a puzzling error when operating on
	// tile compressed images (dev) .. skipping it seems
	// to work OK, error may be related to the bscale reset
	// below that was causing problems
	if (fits_resize_img(is.output_ffp,((is.integer)?LONG_IMG:FLOAT_IMG),naxis,naxes,&is.fits_status))
	  printerror(is.fits_status);
	
	if (0) {
	  // reset bscale  - removed from control because this was screwing things up (unsure why)
	  if (fits_set_bscale(is.output_ffp,1.0,0.0,&is.fits_status))
	    printerror(is.fits_status);
	}
	
	// then write the pixels
	if (fits_write_pix(is.output_ffp,TFLOAT,fpixel,(long)naxes[0]*naxes[1],output_im,&is.fits_status))
	  printerror(is.fits_status);
	// done with this median image.
	if (is.Xray) {
	  // use the just-constructed median and subtract from the
	  // bias subtracted frames, for this extension, for all input files.

	  // open an output file, using CHANNEL or XTENSION
	  // keyword values.
	  //	  file *xfp=fopen("","w");
	  //
	  //	  if (xfp==NULL) {
	  //	    fprintf(stderr,"can't open evlist file.\n");
	  //	    exit(1);
	  //	  }

	  static int idumi=-1;
	  for (ff_ix=0;ff_ix<is.n_ff;ff_ix++) {
	    float *this_imbuf=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
	    if (this_imbuf==NULL) {
	      fprintf(stderr,"can't allocate. (0)\n");
	      usage();
	    }
	    {
	      int hdu_save;
	      fits_get_hdu_num(is.ffp[ff_ix],&hdu_save);
	      if (fits_movabs_hdu(is.ffp[ff_ix],is.edm.data_ext_map[hdu_save],
				  NULL,&is.fits_status) || 
		  fits_read_pix(is.ffp[ff_ix],TFLOAT,
				fpixel,(long)naxes[0]*naxes[1],
				NULL,this_imbuf,NULL,
				&is.fits_status) ||
		  fits_movabs_hdu(is.ffp[ff_ix],hdu_save,
				  NULL,&is.fits_status)) 
		printerror(is.fits_status);
	    }
	    // subtract off this bias
	    m=naxes[0]*naxes[1];
	    while (m--) {
	      this_imbuf[m] -= (is.ibp[ff_ix].bias+output_im[m]);
	    }
	    m=naxes[0]*naxes[1];
	    while (m--) {
	      long addr_ix[]={m%naxes[0],m/naxes[0]};
	      if ((addr_ix[0]-1)*((naxes[0]-2)-addr_ix[0])<0) continue;
	      if ((addr_ix[1]-1)*((naxes[1]-2)-addr_ix[1])<0) continue;
	      // check against event threshold
	      if (this_imbuf[m]<60) continue;
	      float evt[]={this_imbuf[m-naxes[0]-1],
			   this_imbuf[m-naxes[0]],
			   this_imbuf[m-naxes[0]+1],
			   this_imbuf[m-1],
			   this_imbuf[m],
			   this_imbuf[m+1],
			   this_imbuf[m+naxes[0]-1],
			   this_imbuf[m+naxes[0]],
			   this_imbuf[m+naxes[0]+1]};
	      if ((evt[4]> evt[0]) && (evt[4]> evt[1]) && (evt[4]> evt[2]) &&
		  (evt[4]> evt[3]) && 
		  (evt[4]>=evt[5]) && (evt[4]>=evt[6]) && (evt[4]>=evt[7]) &&
		  (evt[4]>=evt[8])) {
		if (is.bias_ffname == NULL) {  // only correct zero values if bias is computed from the stack
		  int j=9;
		  while (j--) {
		    if (evt[j]==0) {
		      evt[j] += 3*gasdev(&idumi);
		    }
		  }
		}
		if (1) {
		  struct data_str ds;
		  ds.mode=0;
		  ds.framenum=ff_ix;
		  ds.chipnum=hdu_ix;
		  ds.x=addr_ix[0];
		  ds.y=addr_ix[1];
		  // if on-the-fly median subtraction is being done *and* the pixel value is zero,
		  // add in some random noise to keep sharp features from being produced in low signal
		  // histograms
		  int j=9;
		  while (j--) {
		    ds.data[j]=floor(evt[j]+0.5);
		  }
		  fwrite(&ds,datastr_size,1,is.xray_output_file);
		} else {
		  fprintf(stderr,"evt at (x,y)=(%ld,%ld):\n",addr_ix[0],addr_ix[1]);
		  fprintf(stderr,"\ttotal = %g\n",floor(evt[0]+evt[1]+evt[2]+
							evt[3]+evt[4]+evt[5]+
							evt[6]+evt[7]+evt[8]));
		  fprintf(stderr,"\t%8g %8g %8g\n",floor(evt[6]),floor(evt[7]),floor(evt[8]));
		  fprintf(stderr,"\t%8g %8g %8g\n",floor(evt[3]),floor(evt[4]),floor(evt[5]));
		  fprintf(stderr,"\t%8g %8g %8g\n",floor(evt[0]),floor(evt[1]),floor(evt[2]));
		}
		// remove from this_imbuf to avoid double counting.
		this_imbuf[m-naxes[0]-1]=0;
		this_imbuf[m-naxes[0]]=0;
		this_imbuf[m-naxes[0]+1]=0;
		this_imbuf[m-1]=0;
		this_imbuf[m]=0;
		this_imbuf[m+1]=0;
		this_imbuf[m+naxes[0]-1]=0;
		this_imbuf[m+naxes[0]]=0;
		this_imbuf[m+naxes[0]+1]=0;
	      }
	    }
	    /*   subtract file/hdu specific bias */
	    /*   and subtract also the bias subtracted median */
	    /*   and look for X-ray events */
	    /* 		      then ship them out. */

	    free(this_imbuf);
	  }
	}
	// append histograms if appropriate.
	if (is.histogram_output) {
	  // go over the image to generate bias and image region histograms and
	  // output them to File *is.histgram_output
	  // last chance before images are freed
	  long min_bias=1000000;
	  long max_bias=-1000000;
	  long min_data=1000000;
	  long max_data=-1000000;
	  long n_bias=0;
	  long n_image=0;
	  long n_exclude=0;
	  
	  long m=naxes[0]*naxes[1];
	  while (m--) {
	    if (is.incl[m] == _bias) {
	      if (min_bias>output_im[m]) min_bias=floor(output_im[m]+0.5);
	      if (max_bias<output_im[m]) max_bias=floor(output_im[m]+0.5);
	      n_bias++;
	    }
	    if (is.incl[m] == _image) {
	      if (min_data>output_im[m]) min_data=floor(output_im[m]+0.5);
	      if (max_data<output_im[m]) max_data=floor(output_im[m]+0.5);
	      n_image++;
	    }
	    if (is.incl[m] == _exclude) n_exclude++;
	  }
	  if (max_data-min_data>2e6) max_data=min_data+2e6;
	  if (max_bias-min_bias>2e6) max_bias=min_bias+2e6;

	  // allocate arrays
	  long *datahist=(long*)malloc((max_data-min_data+1)*sizeof(long));
	  long *biashist=(long*)malloc((max_bias-min_bias+1)*sizeof(long));
	  if ((datahist==NULL) || (biashist==NULL)) {
	    fprintf(stderr,"can't allocate. (10)\n");
	    usage();
	  }
	  // zero out the arrays
	  m=max_bias-min_bias+1;
	  while (m--) biashist[m]=0L;
	  m=max_data-min_data+1;
	  while (m--) datahist[m]=0L;
	  m=naxes[0]*naxes[1];
	  while (m--) {
	    if (is.incl[m] == _bias) {
	      if (floor(output_im[m]+0.5)-min_bias<2e6)
		biashist[(long)(floor(output_im[m]+0.5)-min_bias)]++;
	    }
	    if (is.incl[m] == _image) {
	      if (floor(output_im[m]+0.5)-min_data<2e6)
		datahist[(long)(floor(output_im[m]+0.5)-min_data)]++;
	    }
	  }

	  m=max_bias-min_bias+1;
	  while (m--) {
	    fprintf(is.histogram_output,
		    "%ld %ld\n",m+min_bias,biashist[m]);
	  }
	  fprintf(is.histogram_output,"no no\n");

	  m=max_data-min_data+1;
	  while (m--) {
	    fprintf(is.histogram_output,
		    "%ld %ld\n",m+min_data,datahist[m]);
	  }
	  fprintf(is.histogram_output,"no no\n");
	  fflush(is.histogram_output);
	  free(datahist);
	  free(biashist);
	}
	// free up
	free(im);
	free(output_im);
      }
    }
    free(shipped_out_hdu);
    shipped_out_hdu=NULL;

    if (is.multiplex_mosaic) {
      // combine the already shipped out channels into a MOSAIC extension at the end of the output file.
      // see if there is already a MOSAIC extension and resize, if necessary.
      int num_output_hdus;
      if (fits_get_num_hdus(is.output_ffp,&num_output_hdus,&is.fits_status))
	printerror(is.fits_status);
      int hdu_ix;
      int hdu_type;
      int mosaic_hdu=0;
      int naxis=2;
      long naxes[]={100,100};
      int bitpix=(is.integer)?LONG_IMG:FLOAT_IMG;

      for (hdu_ix=2;hdu_ix<=num_output_hdus;hdu_ix++) { // skip over the PRIMARY hdu
	if (fits_movabs_hdu(is.output_ffp,hdu_ix,&hdu_type,&is.fits_status))
	  printerror(is.fits_status);
	if (hdu_type==IMAGE_HDU) {
	  char extname[2048];
	  if (fits_read_key(is.output_ffp,TSTRING,"EXTNAME",extname,NULL,&is.fits_status)) {
	    if (is.fits_status==KEY_NO_EXIST) {
	      is.fits_status=0;
	    } else {
	      printerror(is.fits_status);
	    }
	  } else {
	    if (strcmp(extname,"MOSAIC")==0) 
	      mosaic_hdu=hdu_ix;
	  }
	} 
      }
      if (mosaic_hdu==0) {
	fprintf(stderr,"didn't find MOSAIC IMAGE_HDU.\ngenerating..\n");
	// existing MOSAIC not found. so create one.
	if (fits_create_img(is.output_ffp,bitpix,naxis,naxes,&is.fits_status))
	  printerror(is.fits_status);
	if (fits_write_key_str(is.output_ffp,"EXTNAME","MOSAIC",NULL,&is.fits_status))
	  printerror(is.fits_status);
	fits_get_hdu_num(is.output_ffp,&hdu_ix);
	mosaic_hdu=hdu_ix;
	if (fits_get_num_hdus(is.output_ffp,&num_output_hdus,&is.fits_status))
	  printerror(is.fits_status);
      }
      // find the expected output mosaic image size from the DETSIZE keyword(s).
      char detsize[4096];
      int found=0;
      for (hdu_ix=2;hdu_ix<=num_output_hdus;hdu_ix++) { 
	// CRAP!! in some files, the entry in hdu=1 contains INCORRECT INFORMATION [1:4352,1:4096]
	// that's pretty lame
	if (found) continue;
	if (hdu_ix==mosaic_hdu) continue;
	if (fits_movabs_hdu(is.output_ffp,hdu_ix,&hdu_type,&is.fits_status))
	  printerror(is.fits_status);
	if (hdu_type!=IMAGE_HDU) continue;
	// look for DETSIZE
	if (fits_read_key(is.output_ffp,TSTRING,"DETSIZE",detsize,NULL,&is.fits_status)) {
	  if (is.fits_status==KEY_NO_EXIST) {
	    is.fits_status=0;
	  } else {
	    printerror(is.fits_status);
	  }
	} else {
	  found=1;
	}
      }
      if (found==0) {
	// detsize keyword doesn't exist in the extensions. last ditch effort to look for it
	// in the primary fits header
	hdu_ix=1;
	if (fits_movabs_hdu(is.output_ffp,hdu_ix,&hdu_type,&is.fits_status))
	  printerror(is.fits_status);
	// look for DETSIZE
	if (fits_read_key(is.output_ffp,TSTRING,"DETSIZE",detsize,NULL,&is.fits_status)) {
	  if (is.fits_status==KEY_NO_EXIST) {
	    is.fits_status=0;
	  } else {
	    if (is.fits_status)
	      printerror(is.fits_status);
	  }
	} else {
	  found=1;
	}
      }
      if (!found) {
	fprintf(stderr,"can't find DETSIZE keyword in input file(s).\n");
	usage();
      }
      int xlim[2],ylim[2];
      sscanf(detsize,"[%d:%d,%d:%d]",&xlim[0],&xlim[1],&ylim[0],&ylim[1]);
      if (fits_movabs_hdu(is.output_ffp,mosaic_hdu,&hdu_type,&is.fits_status))
	printerror(is.fits_status);
      // and resize
      naxes[0]=(xlim[1]-xlim[0]+1);	naxes[1]=(ylim[1]-ylim[0]+1);

      if (fits_resize_img(is.output_ffp,bitpix,naxis,naxes,&is.fits_status))
	printerror(is.fits_status);

      // allocate the mosaic;

      float *mosaic_im=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
      if (mosaic_im==NULL) {
	fprintf(stderr,"can't allocate. (1)\n");
	usage();
      }

      // now go through each segment to construct the mosaic
      int seg_ix;
      for (seg_ix=2;seg_ix<=num_output_hdus;seg_ix++) {
	if (seg_ix==mosaic_hdu) continue;
	if (fits_movabs_hdu(is.output_ffp,seg_ix,&hdu_type,&is.fits_status))
	  printerror(is.fits_status);
	if (hdu_type!=IMAGE_HDU) continue;

	char datasec[2048];
	char detsec[2048];
	int datasec_found=0;
	int detsec_found=0;

	if (fits_read_key(is.output_ffp,TSTRING,"DATASEC",datasec,NULL,&is.fits_status)) {
	  if (is.fits_status==KEY_NO_EXIST) {
	    is.fits_status=0;
	  } else {
	    printerror(is.fits_status);
	  }
	} else {
	  datasec_found=1;
	}

	if (fits_read_key(is.output_ffp,TSTRING,"DETSEC",detsec,NULL,&is.fits_status)) {
	  if (is.fits_status==KEY_NO_EXIST) {
	    is.fits_status=0;
	  } else {
	    printerror(is.fits_status);
	  }
	} else {
	  detsec_found=1;
	}

	if (datasec_found && detsec_found) {
	  int detsec_xlim[2],detsec_ylim[2];
	  int datasec_xlim[2],datasec_ylim[2];
	  sscanf(detsec,"[%d:%d,%d:%d]",
		 &detsec_xlim[0],&detsec_xlim[1],
		 &detsec_ylim[0],&detsec_ylim[1]);
	  sscanf(datasec,"[%d:%d,%d:%d]",
		 &datasec_xlim[0],&datasec_xlim[1],
		 &datasec_ylim[0],&datasec_ylim[1]);
	  // compare segment sizes.
	  if ((abs(detsec_xlim[1]-detsec_xlim[0]) != 
	       abs(datasec_xlim[1]-datasec_xlim[0])) ||
	      (abs(detsec_ylim[1]-detsec_ylim[0]) != 
	       abs(datasec_ylim[1]-datasec_ylim[0]))) {
	    fprintf(stderr,"unexpected/poor correspondence between detsec & datasec:\n");
	    fprintf(stderr,"\tdetsec: %s\n",detsec);
	    fprintf(stderr,"\tdatasec: %s\n",datasec);
	    fprintf(stderr,"attempting a fix.\n");
	    if (detsec_xlim[0]==0) detsec_xlim[0]=1;
	    if (detsec_xlim[1]==0) detsec_xlim[1]=1;
	    if (detsec_ylim[0]==0) detsec_ylim[0]=1;
	    if (detsec_ylim[1]==0) detsec_ylim[1]=1;
	    //usage();
	  }
	  //	  fprintf(stderr,"detsec %s datasec %s\n",detsec,datasec);
	  // datasec is ordered, detsec maybe reversed along either or both axes
	  if (datasec_xlim[1]<datasec_xlim[0]) {
	    int tmp;
	    tmp             = datasec_xlim[1];
	    datasec_xlim[1] = datasec_xlim[0];
	    datasec_xlim[0] = tmp;
	    tmp             = detsec_xlim[1];
	    detsec_xlim[1]  = detsec_xlim[0];
	    detsec_xlim[0]  = tmp;
	    
	  }
	  if (datasec_ylim[1]<datasec_ylim[0]) {
	    int tmp;
	    tmp             = datasec_ylim[1];
	    datasec_ylim[1] = datasec_ylim[0];
	    datasec_ylim[0] = tmp;
	    tmp             = detsec_ylim[1];
	    detsec_ylim[1]  = detsec_ylim[0];
	    detsec_ylim[0]  = tmp;
	  }
	  
	  // allocate & read in segment
	  {
	    int segment_maxdim=5;
	    int segment_bitpix,segment_naxis;
	    long segment_naxes[segment_maxdim];
	    long fpixel[]={1,1};

	    if (fits_get_img_param(is.output_ffp,segment_maxdim,
				   &segment_bitpix,&segment_naxis,segment_naxes,&is.fits_status))
	      printerror(is.fits_status);
	    float *segment_im=(float*)malloc(segment_naxes[0]*segment_naxes[1]*sizeof(float));

	    if (segment_im==NULL) {
	      fprintf(stderr,"can't allocate. (2)\n");
	      usage();
	    }

	    // pixel mapping here shouldn't be done - data is being taken from output.
	    if (fits_read_pix(is.output_ffp,TFLOAT,fpixel,segment_naxes[0]*segment_naxes[1],
			      NULL,segment_im,NULL,&is.fits_status))
	      printerror(is.fits_status);
	    // transfer the datasec onto the detsec for the mosaic.
	    int y_ix,x_ix;
	    int detsec_x_dir=(detsec_xlim[1]>detsec_xlim[0])?1:-1;
	    int detsec_y_dir=(detsec_ylim[1]>detsec_ylim[0])?1:-1;
	    for (y_ix=0;y_ix<(datasec_ylim[1]-datasec_ylim[0]+1);y_ix++) {
	      for (x_ix=0;x_ix<(datasec_xlim[1]-datasec_xlim[0]+1);x_ix++) {
		mosaic_im[(y_ix*detsec_y_dir+detsec_ylim[0]-1)*naxes[0]+
			  (x_ix*detsec_x_dir+detsec_xlim[0]-1)]=
		  segment_im[(y_ix+datasec_ylim[0]-1)*segment_naxes[0]+(x_ix+datasec_xlim[0]-1)];
	      }
	    }
	    free(segment_im);
	  }
	}
      }
      long fpixel[]={1,1};
      if (fits_movabs_hdu(is.output_ffp,mosaic_hdu,&hdu_type,&is.fits_status))
	printerror(is.fits_status);
      if (fits_write_pix(is.output_ffp,TFLOAT,
			 fpixel,(long)naxes[0]*naxes[1],mosaic_im,&is.fits_status))
	printerror(is.fits_status);
      
      free(mosaic_im);
    }
    if (fits_close_file(is.output_ffp,&is.fits_status))
      printerror(is.fits_status);
  }


  free(is.ffp);
  free(is.ffname);
  free(is.hdu);
  free(is.ibp);

  exit(0);
}

int float_compare(const void *a, const void *b) {
  float diff=(*(float*)a)-(*(float*)b);
  if (diff<0) return(-1);
  if (diff>0) return(1);
  return(0);
}

//void compute_image_base_pars(fitsfile *ff,int minser,img_basepars *ibp) {
void compute_image_base_pars(img_stack *isk,int im_ix) {
  fitsfile *ff=isk->ffp[im_ix];
  int maxdim=5;
  int bitpix,naxis;
  long naxes[maxdim];
  int status=0;
  if (fits_get_img_param(ff,maxdim,&bitpix,&naxis,naxes,&status))
    printerror(status);
  // look for DATASEC keyword
  char datasec[2048];
  if (fits_read_key(ff,TSTRING,"DATASEC",datasec,NULL,&status))
    printerror(status);
  //  fprintf(stderr,"got datasec %s\n",datasec);
  int xlim[2],ylim[2];
  sscanf(datasec,"[%d:%d,%d:%d]",&xlim[0],&xlim[1],&ylim[0],&ylim[1]);
  // read in entire image segment, save only bias samples and compute.
  {
    float *img=NULL;
    img=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
    char *incl=isk->incl;
    if (isk->incl == NULL) {
      isk->incl=(char*)malloc(naxes[0]*naxes[1]*sizeof(char));
      incl=isk->incl;
    }

    // arrays for modeling the bias.
    float *serdep_bias=NULL;
    int   *serdep_bias_npix=NULL;
    float *pardep_bias=NULL;
    int   *pardep_bias_npix=NULL;

    serdep_bias=(float*)malloc(naxes[0]*sizeof(float));
    serdep_bias_npix=(int*)malloc(naxes[0]*sizeof(int));
    pardep_bias=(float*)malloc(naxes[1]*sizeof(float));
    pardep_bias_npix=(int*)malloc(naxes[1]*sizeof(int));

    isk->ibp[im_ix].model_bias_array[0].b_lvl = serdep_bias;
    isk->ibp[im_ix].model_bias_array[0].b_nelem = naxes[0];
    isk->ibp[im_ix].model_bias_array[1].b_lvl = pardep_bias;
    isk->ibp[im_ix].model_bias_array[1].b_nelem = naxes[1];
      
    if ((img==NULL) || (incl==NULL) ||
	(serdep_bias==NULL) || (serdep_bias_npix==NULL) ||
	(pardep_bias==NULL) || (pardep_bias_npix==NULL)) {
      fprintf(stderr,"can't allocate. (3)\n");
      usage();
    }

    // first compute the isk->incl mask
    isk->ser_oscan_min = xlim[1] + 4;
    isk->par_oscan_min = ylim[1] + 4;

    isk->ser_oscan_max = naxes[0];
    isk->par_oscan_max = naxes[1];

    isk->bias_model_smooth_npix = 0;

    {
      long m=naxes[0]*naxes[1];
      while (m--) {
	long addr_ix[]={m%naxes[0]+1,m/naxes[0]+1};
	if (((addr_ix[0]-xlim[0])*(xlim[1]-addr_ix[0])>=0) &&
	    ((addr_ix[1]-ylim[0])*(ylim[1]-addr_ix[1])>=0)) {
	  // within DATASEC
	  if (((addr_ix[1]-isk->data_min_par)*
	       (isk->data_max_par-addr_ix[1])>=0) &&
	      ((addr_ix[0]-isk->data_min_ser)*
	       (isk->data_max_ser-addr_ix[0])>=0)) {
	    incl[m] = _image;
	  } else {
	    incl[m] = _exclude;
	  }
	} else {
	  // not within DATASEC
	  if ((addr_ix[0] <= isk->bias_min_ser) ||
	      (addr_ix[1] <= isk->bias_min_par)) {
	    incl[m] = _exclude;
	  } else {
	    incl[m] = _bias;
	  }
	}
      }
    }
    // done with preparation of isk->incl
    
    // zero the serdep & pardep arrays
    {
      int i;
      i=naxes[0]; while (i--) serdep_bias[i]=0.0;
      i=naxes[0]; while (i--) serdep_bias_npix[i]=0;
      i=naxes[1]; while (i--) pardep_bias[i]=0.0;
      i=naxes[1]; while (i--) pardep_bias_npix[i]=0;
    }
	 
    long fpixel[]={1,1};
    {
      int hdu_save;
      fits_get_hdu_num(ff,&hdu_save);
      if (fits_movabs_hdu(ff,isk->edm.data_ext_map[hdu_save],
			  NULL,&status) ||
	  fits_read_pix(ff,TFLOAT,fpixel,naxes[0]*naxes[1],
			NULL,img,NULL,&status) || 
	  fits_movabs_hdu(ff,hdu_save,NULL,&status)) 
	printerror(status);
    }
    long m;
    long min_bias=1000000;
    long max_bias=-1000000;
    long min_data=1000000;
    long max_data=-1000000;
    long latch_max_bias;
    m=naxes[0]*naxes[1];

    
    while (m--) {
      if (incl[m] == _image) {
	if (min_data>img[m]) min_data=floor(img[m]+0.5);
	if (max_data<img[m]) max_data=floor(img[m]+0.5);
      } else {
	long addr_ix[]={m%naxes[0]+1,m/naxes[0]+1};
	// handle the bias modeling histograms slightly differently:
	// check ranges in parallel address for collecting serial dependence
	// check ranges in serial address for collecting parallel dependence
	if ((isk->par_oscan_max - addr_ix[1]) *
	    (addr_ix[1] - isk->par_oscan_min) >= 0) {
	  int j;
	  for (j=addr_ix[0]-1-isk->bias_model_smooth_npix;
	       j<=addr_ix[0]-1+isk->bias_model_smooth_npix;
	       j++) {
	    if ((naxes[0]-1-j)*j>=0) {
	      serdep_bias[j] += img[m];
	      serdep_bias_npix[j] += 1;
	    }
	  }
	}
	if ((isk->ser_oscan_max - addr_ix[0]) *
	    (addr_ix[0] - isk->ser_oscan_min) >= 0) {
	  int j;
	  for (j=addr_ix[1]-1-isk->bias_model_smooth_npix;
	       j<=addr_ix[1]-1+isk->bias_model_smooth_npix;
	       j++) {
	    if ((naxes[1]-1-j)*j>=0) {
	      pardep_bias[j] += img[m];
	      pardep_bias_npix[j] += 1;
	    }
	  }
	}
	if (incl[m] == _bias) {
	  if (min_bias>img[m]) min_bias=floor(img[m]+0.5);
	  if (max_bias<img[m]) {
	    max_bias=floor(img[m]+0.5);
	    latch_max_bias=m;
	  }
	}
      }
    }


    // finish computing the serdep_bias & pardep_bias arrays
    {
      int i;
      int j;
      float bias=0.0;
      int bias_npix=0;

      // first compute mean of the overlap between serdep_bias & pardep_bias
      for (i=xlim[1];i<naxes[0];i++) {
	if (serdep_bias_npix[i]!=0) {
	  bias += serdep_bias[i];
	  bias_npix += serdep_bias_npix[i];
	}
      }
      for (j=ylim[1];j<naxes[1];j++) {
	if (pardep_bias_npix[j]!=0) {
	  bias += pardep_bias[j];
	  bias_npix += pardep_bias_npix[j];
	}
      }
      isk->ibp[im_ix].model_bias = bias/bias_npix;
    }

    // normalize the serdep_bias & pardep_bias arrays since they are refered to by *isk
    {
      int i,j;
      for (i=0;i<naxes[0];i++) {
	serdep_bias[i] /= serdep_bias_npix[i];
	serdep_bias[i] -= isk->ibp[im_ix].model_bias;
      }
      for (j=0;j<naxes[1];j++) {
	pardep_bias[j] /= pardep_bias_npix[j];
	pardep_bias[j] -= isk->ibp[im_ix].model_bias;
      }
    }
    // done with model_bias objects.
    
    //    fprintf(stderr,"min,max bias = (%ld,%ld)\n",min_bias,max_bias);
    //    fprintf(stderr,"min,max data = (%ld,%ld)\n",min_data,max_data);
    // make up histogram

    long *biashist;
    biashist=(long*)malloc((max_bias-min_bias+1)*sizeof(long));
    if (biashist==NULL) {
      fprintf(stderr,"can't allocate. (4)\n");
      usage();
    }
    m=max_bias-min_bias+1;
    while (m--) biashist[m]=0;
    m=naxes[0]*naxes[1];
    while (m--) {
      if (incl[m]==_bias) {
	biashist[(long)(floor(img[m]+0.5)-min_bias)]++;
      }
    }

    // compute bias.. but first, throw out 5 sigma + outliers
    double s_bias=0.0;
    double s2_bias=0.0;
    long n_bias=0;
    m=max_bias-min_bias+1;
    while (m--) {
      n_bias += biashist[m];
      s_bias += biashist[m]*(m+min_bias);
      s2_bias += biashist[m]*pow((m+min_bias),2);
    }
    double tmp_mean = s_bias/n_bias;
    double tmp_rms = sqrt(s2_bias/n_bias - pow(tmp_mean,2));

    s_bias=0.0;
    n_bias=0;
    m=max_bias-min_bias+1;
    while (m--) {
      if (((m+min_bias)-(tmp_mean-5*tmp_rms))*
	  ((tmp_mean+5*tmp_rms)-(m+min_bias))>=0) {
	n_bias += biashist[m];
	s_bias += biashist[m]*(m+min_bias);
      }
    }
    
    isk->ibp[im_ix].bias = s_bias / n_bias;

    if (isk->diagnostic_output) {
      fprintf(isk->diagnostic_output,
	      "! bias dist for image %s, HDU %d bias %g\n",
	      isk->ffname[im_ix],isk->hdu[im_ix],isk->ibp[im_ix].bias);
      fprintf(isk->diagnostic_output,"! max bias loc ix=%ld (ser,par)=(%ld,%ld)\n",
	      latch_max_bias,
	      latch_max_bias%naxes[0]+1,latch_max_bias/naxes[0]+1);
      
      m=max_bias-min_bias+1;
      while (m--) {
	if (((m+min_bias)-(tmp_mean-5*tmp_rms))*
	    ((tmp_mean+5*tmp_rms)-(m+min_bias))>=0) 
	  fprintf(isk->diagnostic_output,"%ld %ld\n",m+min_bias,biashist[m]);
      }
      fprintf(isk->diagnostic_output,"no no\n");
      fflush(isk->diagnostic_output);
    }
    
    free(biashist);

    // now compute the level of this frame after applying a bias correction
    // run through first time to compute offset and then get statistics for 

    // image section
    long *datahist;
    datahist=(long*)malloc((max_data-min_data+1)*sizeof(long));
    if (datahist==NULL) {
      fprintf(stderr,"can't allocate. (5)\n");
      usage();
    }
    m=(max_data-min_data+1);
    while (m--) datahist[m]=0;
    m=naxes[0]*naxes[1];
    while (m--) {
      if (incl[m]==_image) {
	datahist[(long)(floor(img[m]+0.5)-min_data)]++;
      }
    }
    m=max_data-min_data+1;
    if (0) { // compute by mean
      float s_data=0;
      long n_data=0;
      while (m--) {
	// compute the mean to judge gain
	n_data += datahist[m];
	s_data += datahist[m]*(m+min_data);
      }
      isk->ibp[im_ix].gain = (s_data/n_data) - isk->ibp[im_ix].bias;
    } else {
      // compute by median or some sort of truncated mean
      // consider formula that maps between mean & median:
      // lop off outliers, then compute mean of surviving bins
      long n_data=0;
      // find the 50% level as placeholder
      while (m--) {
	n_data += datahist[m];
      }
      // n_data is the number of pixels counted. divide by 2 and count down again
      n_data /= 2;
      m=max_data-min_data+1; // start at the top again
      while (m-- && (n_data>0)) {
	n_data -= datahist[m];
      }
      // use m+min_data as the median
      isk->ibp[im_ix].gain = ( m + min_data ) - isk->ibp[im_ix].bias;
    }
    free(datahist);

    if (isk->img_basepar_output) {
      fprintf(isk->img_basepar_output,"filename %s HDU %d bias %g scale %g\n",
	      isk->ffname[im_ix],
	      isk->hdu[im_ix],
	      isk->ibp[im_ix].bias,
	      isk->ibp[im_ix].gain);
      fflush(isk->img_basepar_output);
    }
    free(img);
  }
}

void get_ibp_list(img_stack *isk) {

  char myline[4096];
  char myfile[2048];
  int  ibp_chunk=4;
  int  ibp_n=0;
  int  *ibp_hdu=(int*)malloc(ibp_chunk*sizeof(int));
  float *ibp_bias=(float*)malloc(ibp_chunk*sizeof(float));
  float *ibp_scale=(float*)malloc(ibp_chunk*sizeof(float));
  if ((ibp_hdu==NULL) || (ibp_bias==NULL) || (ibp_scale==NULL)) {
    fprintf(stderr,"can't allocate. (6)\n");
    usage();
  }
    
  int myhdu;
  float mybias,myscale,mytweak;
  while (fgets(myline,2048,isk->img_basepar_input)) {
    int nret=sscanf(myline,"%*s %s %*s %d %*s %f %*s %f %f\n",
		    myfile,&myhdu,&mybias,&myscale,&mytweak);
    if (nret==5) {
      fprintf(stderr,"tweaking.. myscale %g -> %g\n",myscale,myscale*mytweak);
      myscale *= mytweak;
    }
    if (ibp_n >= ibp_chunk) { //realloc if necessary
      ibp_chunk *= 2;
      ibp_hdu=(int*)realloc(ibp_hdu,ibp_chunk*sizeof(int));
      ibp_bias=(float*)realloc(ibp_bias,ibp_chunk*sizeof(float));
      ibp_scale=(float*)realloc(ibp_scale,ibp_chunk*sizeof(float));
      if ((ibp_hdu==NULL) || (ibp_bias==NULL) || (ibp_scale==NULL)) {
	fprintf(stderr,"can't allocate. (7)\n");
	usage();
      }
    }
    ibp_hdu[ibp_n]=myhdu;
    ibp_bias[ibp_n]=mybias;
    ibp_scale[ibp_n]=myscale;
    ibp_n++;
  }
  fclose(isk->img_basepar_input);
  // on second thought, don't reset this to zero. it can convey information downstream.
  // isk->img_basepar_input=NULL; 
  // identify largest HDU
  int ibp_ix;
  int max_hdu=0;
  ibp_ix=ibp_n;
  while (ibp_ix--) {
    if (max_hdu<ibp_hdu[ibp_ix]) max_hdu=ibp_hdu[ibp_ix];
  }
  float *ibp_bias_by_hdu=(float*)malloc((max_hdu+1)*sizeof(float));
  float *ibp_scale_by_hdu=(float*)malloc((max_hdu+1)*sizeof(float));
  int   *ibp_n_by_hdu=(int*)malloc((max_hdu+1)*sizeof(int));
  if ((ibp_bias_by_hdu==NULL) ||
      (ibp_scale_by_hdu==NULL) ||
      (ibp_n_by_hdu==NULL)) {
    fprintf(stderr,"can't allocate. (8)\n");
    usage();
  }
  int hdu_ix=max_hdu+1;
  while (hdu_ix--) 
    ibp_bias_by_hdu[hdu_ix]=ibp_scale_by_hdu[hdu_ix]=ibp_n_by_hdu[hdu_ix]=0;

  ibp_ix=ibp_n;
  while (ibp_ix--) {
    ibp_bias_by_hdu[ibp_hdu[ibp_ix]]  += ibp_bias[ibp_ix];
    ibp_scale_by_hdu[ibp_hdu[ibp_ix]] += ibp_scale[ibp_ix];
    ibp_n_by_hdu[ibp_hdu[ibp_ix]]++;
  }
  hdu_ix=max_hdu+1;
  while (hdu_ix--) {
    ibp_bias_by_hdu[hdu_ix]  /= ibp_n_by_hdu[hdu_ix];
    ibp_scale_by_hdu[hdu_ix] /= ibp_n_by_hdu[hdu_ix];
  }
  // finally hand over results to the isk structure
  isk->bias_by_hdu=ibp_bias_by_hdu;
  isk->scale_by_hdu=ibp_scale_by_hdu;
  // and free up allocated space that's safe to deallocate.
  free(ibp_n_by_hdu);
  free(ibp_scale);
  free(ibp_bias);
  free(ibp_hdu);
  return;
}

void printerror( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
 
  char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
 
  if (status)
    fprintf(stderr, "\n*** Error occurred during program execution ***\n");
 
  fits_get_errstatus(status, status_str);   /* get the error description */
  fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);
 
  /* get first message; null if stack is empty */
  if ( fits_read_errmsg(errmsg) )
    {
      fprintf(stderr, "\nError message stack:\n");
      fprintf(stderr, " %s\n", errmsg);
 
      while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
        fprintf(stderr, " %s\n", errmsg);
    }
 
  exit( status );       /* terminate the program, returning error status */
  usage();
}

void img_isr(img_stack *isk,int ff_ix,float *im,long *naxes) {
  long fpixel_read[]={1L,1L};
  long npixel_read=naxes[0]*naxes[1];
  img_isr_chunk(isk,ff_ix,im,naxes,fpixel_read,npixel_read);
}

void img_isr_chunk(img_stack *isk,int ff_ix,float *im,long *naxes,
		   long *fpixel_read,long npixel_read) {
  
  long pix_ix_offset = (fpixel_read[1]-1)*naxes[0] + (fpixel_read[0]-1);
  long pix_ix = npixel_read;
  
  while (pix_ix--) {
    if (isk->bias_1d_plus_1d_model) {
      long addr_ix[] = {(pix_ix+pix_ix_offset)%naxes[0],
			(pix_ix+pix_ix_offset)/naxes[0]};
      long par_addr=addr_ix[1];
      if (par_addr==0) {
	im[pix_ix] -= (isk->ibp[ff_ix].model_bias_array[0].b_lvl[addr_ix[0]]
		       + isk->ibp[ff_ix].model_bias_array[1].b_lvl[addr_ix[1]]
		       + isk->ibp[ff_ix].model_bias );
      } else {
	float u = (addr_ix[0])/(1.0*naxes[0]); // weight for bias at par_addr. (1-u) for par_addr-1.
	im[pix_ix] -= (isk->ibp[ff_ix].model_bias_array[0].b_lvl[addr_ix[0]]
		       + (u*isk->ibp[ff_ix].model_bias_array[1].b_lvl[addr_ix[1]]
			  + (1-u)*isk->ibp[ff_ix].model_bias_array[1].b_lvl[addr_ix[1]-1])
		       + isk->ibp[ff_ix].model_bias );
      }
    } else { // must be is.bias_subtract
      im[pix_ix] -= isk->ibp[ff_ix].bias;
    }
    if (isk->normalize)
      im[pix_ix] /= isk->ibp[ff_ix].gain;
  }
}

void usage (void) {
  char *usg[]=
    {"usage: mef_combine [options] (required input) mef_0 mef_1 mef_2 .. mef_n",
     " required input: ",
     "   (-m, --median | -a, --average | -s, --sum | -D, --diff | -X, --Xray)",
     " options: ",
     "   -o, --output_file\t<filename> [mef_combine probably won't do anything",
     "\t\t\tuseful without specifying an output (mef) file.]",
     "   -v, --verbose [not much helpful information now]",
     "   -z, --normalize [output will be scaled (useful for multiplexed images)]",
     "   -i, --integer [will store output as a long image instead of floating point]",
     "   -M, --multiplex [appends output file with a multiplexed (mosaic) image]",
     "   -b, --bias_subtract [computed & subtracted on a frame by frame basis]",
     "   --bias_file\t<bias_filename> [specify file to use as bias]",
     "   --88 [remaps improperly ordered/labeled data seen in CCD250-088]",
     "   --bias_min_parallel\t<mimimum parallel address for bias sampling>",
     "   --bias_min_serial\t<miminum serial address for bias sampling>",
     "   --data_min_parallel\t<minimum parallel address for data sampling>",
     "   --data_max_parallel\t<maximum parallel address for data sampling>",
     "   --data_min_serial\t<minimum serial address for data sampling>",
     "   --data_max_serial\t<maximum serial address for data sampling>",
     "   --diagnostic_output\t<filename> [text file of bias & data distributions]",
     "   --img_basepar_output\t<filename> [can be reused as input file below",
     "\t\t\tand can speed up subsequent execution. ",
     "\t\t\tIf --img_basepar_input is specified at the same time",
     "\t\t\tthis output file will reflect a unique set of",
     "\t\t\tbasepars which can be appended with \"tweaks\"",
     "\t\t\tfor subsequent execution]",
     "   --img_basepar_input\t<filename> [file containing list of biases and gains",
     "\t\t\tthat will be used in place of on-the-fly estimation.",
     "\t\t\tUsually this would be generated using the ",
     "\t\t\t--img_basepar_output option above.]",
     "   --Xray_output_file_prefix <prefix> [filename prefix for X-ray (evlist) data]"
    };
  int usage_lines=sizeof(usg)/sizeof(char*);
  int i=-1;
  while(++i<usage_lines) {
    fprintf(stderr,"%s\n",usg[i]);
  }
  fprintf(stderr,"\n");
  exit(1);
}
