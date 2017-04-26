#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>

int main()
{
	int status;
	int pid = fork();
	int i=0;
	char *filename = malloc(sizeof(char*)*100);

	if(pid==0)  //producer
	{
		while(i<100){		
		sprintf(filename,"test_%ld.csv",time(NULL));
		printf("\tWriting %d.\n",i);
		FILE* fp = fopen(filename,"w");
		printf("Opened file %s\n",filename);
		fprintf(fp,"%d",i);
		fclose(fp);
		i++;
		usleep(1000000);
		}
	}
	else if(pid>0)  //consumer
	{
		while(1){
		system("ls test_*.csv > filename_list.txt");
		printf("Generated filename_list.txt\n");
		usleep(500000);
		/* reading
		FILE* fp = fopen("test.txt","r");
		if(fp==NULL)
			printf("Failed opening\n");
		else{
			printf("Sucessfully opened.\n");
			char *buf = malloc(sizeof(char*)*100);
			fscanf(fp, "%s", buf);
			fprintf(stdout,"\tRead file content:%s\n",buf);
			usleep(500000);
		}
		fclose(fp);
		*/
		}		
	}
	else
		printf("Fork failed.\n");

	do
	{
		wait(&status);
		printf("waiting..\n");
	}while(status==-1);// <<<<<--wait for child
	exit(0);

}
