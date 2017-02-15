all: watermark 

watermark: watermark.cpp
	g++ watermark.cpp -o watermark -lpng -I.

clean:
	rm watermark 
