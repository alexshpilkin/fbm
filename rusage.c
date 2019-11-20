#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <sysexits.h>
#include <unistd.h>

int main(int argc, char **argv) {
	pid_t pid;
	int status;
	struct rusage ru;

	if (argc < 2) {
		fprintf(stderr, "usage: %s COMMAND ...\n", argv[0]);
		return EX_USAGE;
	}
	if ((pid = fork()) < 0) {
		perror(argv[0]);
		return EXIT_FAILURE;
	}
	if (!pid)
		execvp(argv[1], argv + 1);

	if (waitpid(pid, &status, 0) < 0) {
		perror(argv[0]);
		return EXIT_FAILURE;
	}
	getrusage(RUSAGE_CHILDREN, &ru);
	printf("# usage %li.%.6li %li\n",
	       (long)ru.ru_utime.tv_sec,
	       (long)ru.ru_utime.tv_usec,
	       ru.ru_maxrss);
	return WIFEXITED(status) ? WEXITSTATUS(status) : EXIT_FAILURE;
}
