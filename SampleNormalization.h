#ifndef SAMPLE_NORMALIZATION_H
#define SAMPLE_NORMALIZATION_H

namespace dpl{

double **ReadSample( char *FileName, int sampleNumber, int sampleElementNumber ){
	FILE *fp;
	fp = fopen( FileName, "r");
        if( fp == NULL ){
		printf("Could not find sample file.\n");		
		exit(0);
	}

	double **sample = (double**)malloc(sampleNumber*sizeof(double*));
        for( unsigned int i=0; i<sampleNumber; i++ )
		sample[i] = (double*)malloc(sampleElementNumber*sizeof(double));
	
	for( unsigned int i=0; i<sampleElementNumber; i++ ){
		for( unsigned int j=0; j<sampleNumber; j++ )
        		fscanf(fp, "%lf", &sample[j][i]); 
	}
	fclose(fp);
	return sample;
}

void SampleNormalization( double **sample, int sampleNumber, int sampleElementNumber )
{ 	
	for( unsigned int i=0; i<sampleNumber; i++ ){
		
		double mean = 0;
		for( unsigned int j=0; j<sampleElementNumber; j++ )
             		mean += sample[i][j];

		mean = mean/sampleElementNumber;

		double sum = 0;	
		for( unsigned int j=0; j<sampleElementNumber; j++ ){
             		sample[i][j] -= mean;
			sum += sample[i][j]*sample[i][j];
		}
		
		double sqrtSum = sqrt(sum);
		if( sqrtSum!=0 ){
			for( unsigned int j=0; j<sampleElementNumber; j++ )
				sample[i][j] = (sample[i][j]/sqrtSum);
		}
	}	
}

std::vector<std::string> split(std::string const &input) { 
	std::istringstream buffer(input);
	std::vector<std::string> ret;
	std::copy(std::istream_iterator<std::string>(buffer),
	std::istream_iterator<std::string>(),
	std::back_inserter(ret));
	return ret;
}

int getSampleNumber( char *sampleFileName ) {
	std::string line;
  	std::ifstream file (sampleFileName);
	std::vector<std::string> firstLine;
  	if (file.is_open()){
    		if( getline (file,line) )
      			firstLine = split(line);
    		file.close();
  	}
	return firstLine.size();
}

void clearSample( int sampleNumber, double **sample ){
	
	for( unsigned int i=0; i<sampleNumber; i++ )
		free(sample[i]);
	free(sample);
}


}
#endif /* Sample Normalization */
