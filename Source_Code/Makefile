CC = gcc
CFLAG = -lmraa -lfann

default: stride

stride: main.c classifiers.h
	$(CC) $(CFLAGS) -o $@ main.c classifiers.c process.c collect_training_data.c collection.c continuous_data.c process_testing_data.c test_neural_network.c test_sample.c training_file.c

clean:
	rm -f stride_detection
