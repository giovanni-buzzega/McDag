# mcs_power_dag

Compiling: 
```shell
g++ -o mcdag mcdag.cpp
```

usage:

	./mcdag [FLAG...] [STRING...]


Flags:

	-k num_strings		set the number of strings to compute mcs on
	-n len_strings		set the length of the strings
	-s alphabet_size	set the alphabet size
	-m			print all mcs (!!)
	-d			print McDag
	-l			print the distribution of the mcs lengths
	-r			print the number of lcs-1 and lcs
	-z			minimalize McDag
	-a			apply all flags to approximate indices
	-i			do not use McDag optimizations (slow)
	-h, --help		show this help


