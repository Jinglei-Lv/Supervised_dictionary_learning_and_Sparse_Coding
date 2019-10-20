#ifndef DICTIONARY_GENERATION_H
#define DICTIONARY_GENERATION_H

namespace dpl{

double **InitializeDictionary( int featureNumber, int sampleElementNumber ){

	double **Wd = (double**)malloc(sampleElementNumber*sizeof(double*));	
	for( unsigned int i=0; i<sampleElementNumber; i++ )
		Wd[i] = (double*)malloc(featureNumber*sizeof(double));	

	return Wd;	
}

double **GenerateRandomDictionary( int featureNumber, int sampleElementNumber ){

	double **Wd = InitializeDictionary( featureNumber, sampleElementNumber );		
	srand((unsigned)time(0));
	unsigned int myseed = (unsigned int) RAND_MAX * rand();	
	
	for( unsigned int i=0; i<sampleElementNumber; i++ ){
		for( unsigned int j=0; j<featureNumber; j++ ){
			Wd[i][j] = 2*(double) (rand() / (RAND_MAX + 1.0))-1;
		}
	}
	return Wd;
}

double **GenerateRandomPatchDictionary( int featureNumber, int sampleElementNumber, int sampleNumber, double **sample ){

	double **Wd = InitializeDictionary( featureNumber, sampleElementNumber );
	srand((unsigned)time(0));
	unsigned int myseed = (unsigned int) RAND_MAX * rand();

	for( unsigned int i=0; i<featureNumber; i++ ){
        	unsigned int index = rand()%sampleNumber;
        	for( unsigned int j=0; j<sampleElementNumber; j++ )
            		Wd[j][i]=sample[index][j];
	}
	return Wd;
}

void DictionaryNormalization( int featureNumber, int sampleElementNumber, double **Wd ){
	for( unsigned int i=0; i<featureNumber; i++ ){
		double sum = 0;
		for( unsigned int j=0; j<sampleElementNumber; j++ )
			sum += Wd[j][i]*Wd[j][i];
		sum = sqrt(sum);

		if( sum!=0 ){
			for( unsigned int j=0; j<sampleElementNumber; j++ )
				Wd[j][i] = Wd[j][i]/sum;
		}
	}
}

double **readDictionary( char *FileName, int featureNumber, int sampleElementNumber ) {

	double **Wd = InitializeDictionary( featureNumber, sampleElementNumber );
	FILE *fp;
    	fp = fopen( FileName, "r");
    	if( fp == NULL ){
		printf("could not find dictionary file %s\n", FileName);
        	exit(0);
	}
    	for( unsigned int i=0; i<sampleElementNumber; i++ ){
        	for( unsigned int j=0; j<featureNumber; j++)
	        	fscanf(fp, "%lf", &Wd[i][j]);		
	}
	fclose(fp);
	return Wd;
}

void saveDictionary( int featureNumber, int sampleElementNumber, double **Wd, char *dictionaryName ){

	FILE *fp;
        fp = fopen( dictionaryName, "w");
        if( fp == NULL ){
		printf("could not find dictionary file: %s\n", dictionaryName);
            	exit(0);
	}

	for( unsigned int i=0; i<sampleElementNumber; i++ ){
		for( unsigned int j=0; j<featureNumber; j++)
	        	fprintf(fp, "%.50lf ", Wd[i][j]);
		fprintf(fp, "\n");		
	}
	fclose(fp);
	printf("Save dictionary file in %s\n", dictionaryName);
}


void clearDictionary( int sampleElementNumber, double **Wd ){
	for( unsigned int i=0; i<sampleElementNumber; i++ ){
		free(Wd[i]);	
	}	
	free(Wd);
}

}

#endif /* Dictionary Generation*/
