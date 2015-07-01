PAGES=$(wildcard src/*.page)
SRC=$(wildcard ../../include/tut/*.hpp)

all: webgen doxygen

webgen: src/default.css src/default.template src/virtual ${PAGES}
	@webgen

doxygen: Doxyfile dox/header.html dox/footer.html dox/doxygen.css dox/main.dox ${SRC}
	@doxygen

deploy: all
	scp -r out/* khorniszon,tut-framework@web.sourceforge.net:/home/groups/t/tu/tut-framework/htdocs

clean:
	@rm -rf out
	@rm -f doxywarn.txt
	@rm -f webgen.cache
