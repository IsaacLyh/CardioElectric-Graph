all:
	clang++ -std=c++11 -stdlib=libc++ solve.cpp;
	clang++ -std=c++11 -stdlib=libc++ Plotting.cpp;
	clang++ -std=c++11 -stdlib=libc++ cmdLine.cpp;
	clang++ -std=c++11 -stdlib=libc++ Report.cpp;
	clang++ -std=c++11 -stdlib=libc++ utils.cpp;
	clang++ -std=c++11 -stdlib=libc++ helper.cpp;
	clang++ -std=c++11 -stdlib=libc++ Timer.cpp;
	clang++ -std=c++11 -stdlib=libc++ apf.cpp;

clean:
	rm *.o;
	rm apf;
	rm *.core;
