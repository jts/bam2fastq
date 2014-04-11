CXXFLAGS += -I samtools -O3
LDFLAGS += -lbam -Lsamtools -lz -pthread

SRC = bam2fastq.cpp
BAM = samtools/libbam.a
AUX = LICENSE Makefile README.txt HISTORY.txt

OBJ = $(SRC:%.cpp=%.o)

bam2fastq: $(OBJ) $(BAM)
	$(CXX) $(OBJ) $(LDFLAGS) $(CXXFLAGS) -o bam2fastq

.PHONY: clean
clean:
	$(RM) $(OBJ)
	$(RM) bam2fastq

.PHONY: cleanall
cleanall: clean
	$(MAKE) -C samtools clean
    
$(BAM):
	$(MAKE) -C samtools lib

.PHONY: dist
dist: VERSION := `sed -e '/version\[\]/!d' \
	           -e 's/.*"\([^"]*\)".*/\1/' \
			   -e 'q' bam2fastq.cpp`
dist: NAME := bam2fastq-$(VERSION)
dist:
	@echo "Creating scratch directory"
	@mkdir $(NAME)
	@echo "Linking source"
	@cd $(NAME) && ln -s $(SRC:%=../%) $(AUX:%=../%) ../samtools .
	@echo "Creating tarball"
	@tar -czhf $(NAME).tgz --exclude='*.o' --exclude='*.a' $(NAME)
	@echo "Removing scratch directory"
	@rm -rf $(NAME)
