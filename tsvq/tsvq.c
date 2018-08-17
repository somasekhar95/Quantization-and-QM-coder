#include "tsvq.h"

char      *programname;
int       dim;
DISTTYPE  mult_offset;
DISTTYPE  thresh;

main(argc, argv)
     int argc;
     char *argv[];
{
  TreeNode  *root;               /* root of codebook tree */
  SlopeList *(*get_next_node)(); /* pointer to maxnode() and top_of_list() */
  SlopeList *nextslope;          /* the pointer to the node to be replaced */
  FILE      *trsqfile;           /* blocked training sequence file */
  FILE      *codebookfile;       /* file to contain codebook */
  FILE      *statfile;           /* file to contain counts and distortions */
  char      trsqname[NAME_MAX];  /* name of training sequence file */
  char      codebookname[NAME_MAX]; /* name of codebook file */
  char      statname[NAME_MAX];  /* name of stat file */
  char      option;              /* used for command lint interpretation */
  char      *cp;                 /* used for command lint interpretation */
  DISTTYPE  distortion;          /* distortion of codebook and training set */
  DISTTYPE  rate,stoppingrate;   /* current rate, and rate to stop algorithm */
  DISTTYPE  tempdist;            /* temporary variable */
  double    entropy;              /* the empirical entropy of the codebook */
  float     userrate;            /* rate needs to be read as float */
  float     userthresh;          /* a number from 0 to 1 */
  float     useroffset;
  long      bits;                /* number bits used to encode training set */
  long      trsqcount;           /* number of training vectors */
  long      numnodes;            /* number of nodes in the tree */
  BOOLEAN  (*insert_slope)();    /* forced_insert() and conditional_insert() */
  BOOLEAN   balanceflag;         /* if true, the tree is balanced */

  /* assign defaults */
  strcpy(trsqname, DEF_trsqname);
  strcpy(codebookname, DEF_codebookname);
  strcpy(statname, DEF_statname);

  programname = *argv;
  dim = DEF_dim;
  mult_offset = DEF_mult_offset;
  userthresh = DEF_thresh;
  thresh = (DISTTYPE) userthresh;
  userrate = DEF_rate;
  stoppingrate = (DISTTYPE) userrate;
  balanceflag = FALSE;

  /* if no options entered, list all of the defaults */
  if (argc == 1) {
    printf("%s %s %s\n",USAGE,programname,HOWTOUSE_TSVQ);
    printf("\nOPTIONS   DESCRIPTIONS                         DEFAULTS\n");
    printf("-t        blocked training sequence            %s\n",trsqname);
    printf("-c        codebook                             %s\n",codebookname);
    printf("-s        codebook statistics                  %s\n",statname);
    printf("-d        vector dimension                     %d\n",dim);
    printf("-r        stopping rate (bits per vector)      %g\n",userrate);
    printf("-m        multiplicative offset                %g\n",mult_offset);
    printf("-h        convergence threshold                %g\n",userthresh);
    printf("-B        produce a balanced tree\n");
    printf("\n");
    fflush(stdout);
    exit(0);
  }

  /* read and interpret command line arguments */
  while (--argc && ++argv) {
    if (*argv[0]=='-' && strlen(*argv)==2) { /* each option has 1 letter */
      option = *++argv[0];
      if (option == 'B') {balanceflag = TRUE;} /* examine the flag */
      else if (--argc && ++argv) { /* examine the option */
        switch(option) { /* examine the option letter */
        case 't':
          strncpy(trsqname,*argv,NAME_MAX - 10);
          break;
        case 'c':
          strncpy(codebookname,*argv,NAME_MAX);
          break;
        case 's':
          strncpy(statname,*argv,NAME_MAX);
          break;
        case 'd':
	  sscanf(*argv, "%d", &dim);
          break;
        case 'r':
	  sscanf(*argv, "%f", &userrate);
	  stoppingrate = ((DISTTYPE) userrate);
          break;
        case 'm':
	  sscanf(*argv, "%f", &useroffset);
	  mult_offset = (DISTTYPE) useroffset;
          break;
        case 'h':
	  sscanf(*argv, "%f", &userthresh);
	  thresh = ((DISTTYPE) userthresh);
          break;
        default:
          fprintf(stderr,"%s: %c: %s\n",programname,option,NOTOPTION);
          exit(1);
          break;
        }
      }
      else {

        fprintf(stderr,"%s %s %s\n",USAGE,programname,HOWTOUSE_TSVQ);
        exit(2);
      }
    }

    else if (*argv[0] == '-') { /* user entered unknown option */
      ++argv[0];
      fprintf(stderr,"%s: %s: %s\n",programname,*argv,NOTOPTION);
      exit(3);
    }
    else { /* user entered unknown string */
      fprintf(stderr,"%s: %s: %s\n",programname,*argv,NOTOPTION);
      exit(4);
    }
  }

  /* user entered invalid arguments */
  if (dim <= 0 || userrate < 0 || mult_offset <= 0 || thresh < 0 || thresh>1){
    fprintf(stderr,"%s %s %s\n",USAGE,programname,HOWTOUSE_TSVQ);
    exit(5);
  }

  /* user did not enter an file names */
  if (strlen(trsqname) == 0) {
    fprintf(stderr,"%s %s %s\n",USAGE,programname,HOWTOUSE_TSVQ);
    exit(6);
  }
  if (strlen(codebookname) == 0) {
    sprintf(codebookname,"%s%s",trsqname,DEF_APPEND_CDBK);
  }
  if (strlen(statname) == 0) {
    sprintf(statname,"%s%s",trsqname,DEF_APPEND_STAT);
  }

  /* user entered an input name which is the same as an output name */
  if (strncmp(cp = codebookname,trsqname,NAME_MAX) == 0 ||
      strncmp(cp = statname,trsqname,NAME_MAX) == 0) {
    fprintf(stderr,"%s: %s %s %s %s: %s\n",
            programname,cp,AND,trsqname,ARESAME,ABORT_TSVQ);
    exit(7);
  }
  if (strncmp(statname,codebookname,NAME_MAX) == 0) {
    fprintf(stderr,"%s: %s %s %s %s: %s\n",
            programname,statname,AND,trsqname,ARESAME,ABORT_TSVQ);
    exit(8);
  }

  /* open files */
  if(!(trsqfile = fopen(cp = trsqname,"r")) ||
     !(codebookfile = fopen(cp = codebookname,"w")) ||
     !(statfile = fopen(cp = statname,"w"))) {
    fprintf(stderr,"%s: %s: %s\n",programname,cp,NOTFOUND);
    exit(9);
  }

  printf("\nOPTIONS   DESCRIPTIONS                         SETTINGS\n");
  printf("-t        blocked training sequence            %s\n",trsqname);
  printf("-c        codebook                             %s\n",codebookname);
  printf("-s        codebook statistics                  %s\n",statname);
  printf("-d        vector dimension                     %d\n",dim);
  printf("-r        stopping rate (bits per vector)      %g\n",stoppingrate);
  printf("-m        multiplicative offset                %g\n",mult_offset);
  printf("-h        convergence threshold                %g\n",thresh);
  if(balanceflag) {
    printf("-B        produce a balanced tree\n");
  }
  /* assign the appropriate slope node search function*/
  if(balanceflag) {
    get_next_node = find_oldest_entry;
    insert_slope = forced_insert;
  }
  else {
    get_next_node = find_maxslope;
    insert_slope = conditional_insert;
  }

  /* create and initialize the root of the main tree, it there is no data,
     then initialize_tree will return FALSE */
  if(!(root = create_root())) exit(10);
  if(!(initialize_tree(root,trsqfile,trsqname))) exit(11);

  fclose(trsqfile);

  trsqcount = root->count;
  distortion = root->avmse * (DISTTYPE) trsqcount;
  bits = 0;
  rate = ((DISTTYPE) bits) / ((DISTTYPE) trsqcount);
  tempdist = distortion / ((DISTTYPE) trsqcount);

  printf("          number of training vectors           %d\n",root->count);
  printf("\n");
  printf("           RATE          DISTORTION\n");
  printf("%15f     %15f\n", rate, tempdist);
  fflush(stdout);

  /* initialize the slopelist, split the root node and do the GLA
     on the root's children, insert the root into the slopelist */
  if(!(initialize_slopelist())) { exit(12); }
  if(!(lloyd(root))) { exit(13); }
  if(!(insert_slope(root))) { exit(14); }

  /* until the stopping rate is achieved, grow the tree by adding to the
     tree the children of the node with the largest slope. then split
     those children and add the children to the slopelist */
  while (rate < stoppingrate) {

    /* find the node with the largest slope */
    if(!(nextslope = get_next_node())) {
      exit(15);
    }

    /* if the largest slope is 0 the tree cannot be grown anymore */
    if (nextslope->slope <= 0.0) {
      printf("%s: %s\n",programname,HALT_TSVQ);
      fflush(stdout);
      break;
    }

    /* include nextslope node's children in the codebook tree */
    nextslope->node->left_child->designed = TRUE;
    nextslope->node->right_child->designed = TRUE;

    /* print the rate and distortion information */
    update_rate(nextslope,&bits,&distortion);
    rate = ((DISTTYPE) bits) / ((DISTTYPE) trsqcount);
    tempdist = distortion / ((DISTTYPE) trsqcount);
    printf("%15f     %15f\n", rate, tempdist);
    fflush(stdout);

    /* split the nextslope's children if possible and add to the list
       to be considered for splitting if the slope is large enough */

    if(!(lloyd(nextslope->node->left_child))) { exit(16); }
    if(!(insert_slope(nextslope->node->left_child))) { exit(17); }

    if(!(lloyd(nextslope->node->right_child))) { exit(18); }
    if(!(insert_slope(nextslope->node->right_child))) { exit(19); }

    /* remove the nextslope from further consideration */
    delete_slopelist(nextslope);
  }

  /* find the number of nodes and print data to output files */
  if((numnodes = tree_clean(root)) == 0) exit(20);
  if((entropy = empirical_entropy(root,numnodes)) < 0) exit(21);
  printf("\nThe number of nodes is %ld\n",numnodes);
  printf("The empirical entropy is %lf bits per vector\n\n",entropy);
  if(!(write_codebook(root,numnodes,codebookfile,codebookname))) {
    exit(22);
  }
  if(!(write_stat(root,numnodes,statfile,statname))) {
    exit(23);
  }
  exit(0);
}
