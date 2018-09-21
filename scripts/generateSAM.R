library(readr)
library(R6)

options(scipen=999)

FakeSAM <- R6Class("FakeSAM",
    public=list(
        initialize = function(filename) {
            # Open sam file.
            private$filename = filename
            private$con = file(filename, "w")            
            
            # Read header to get chromosome sizes.
            all_header = scan("input/header.txt", character(), sep="\n")[-1]
            all_chr = gsub(".*SN:(\\w+).LN:(.*)", "\\1", all_header)
            private$chr_size = as.integer(gsub(".*SN:(\\w+).LN:(.*)", "\\2", all_header))
            names(private$chr_size) = all_chr
            
            # Write SAM header.
            header_lines <- read_file("input/header.txt")
            cat(header_lines, file=private$con)
        },
        add_flat_coverage = function(seqname, start_position, end_position, coverage_value) {
            read_length = 50
            
            # Add the coverage information to the summary data-frame.
            cov_df = data.frame(seqnames=seqname, 
                                start=start_position, 
                                end=end_position, 
                                raw=coverage_value)
            private$coverage = rbind(private$coverage, cov_df)
                               
            # We only support 50-nt reads for now.
            stopifnot((end_position - start_position) %% read_length == 0)
            
            width_nt = end_position - start_position
            width_reads = width_nt / read_length
            total_reads = width_reads * coverage_value
            
            # Generate random read numbers
            read_ids = ceiling(runif(total_reads, 0, 999999999999))
            
            # Loop until we've generated all reads.
            read_count=1
            for(i in 1:width_reads) {
                for(j in 1:coverage_value) {
                    read_start = start_position+((i-1) * read_length)
                    read_end = read_start + read_length - 1
                    
                    this_align = private$align_line       
                    this_align = gsub("READNAME", paste0("symcoverage.", read_ids[read_count]), this_align)
                    this_align = gsub("POSITION", read_start, this_align)
               
                    cat(this_align, file=private$con)
                    read_gr = GRanges(paste0(seqname, ":", read_start, "-", read_end))
                    private$all_ranges = c(private$all_ranges, read_gr)
                    private$read_count = private$read_count + 1
                }
            }
        },
        write_file_and_close = function() {
            # Close the SAM file.
            close(private$con)
            
            # Calculate RPM coverages for each interval, and write operation summary..
            private$coverage$rpm = private$coverage$raw * 1000000 / private$read_count
            write.table(private$coverage, file=paste0(private$filename, ".coverage.df"))
            
            # Calculate coverages based on read intervals.
            cov_rle = coverage(private$all_ranges)
            
            # Pad coverages to account for chromosome sizes.
            for(chr in names(cov_rle)) {
                padding_zero = Rle(0, private$chr_size[[chr]] - sum(runLength(cov_rle[[chr]])))
                cov_rle[[chr]] = c(cov_rle[[chr]], padding_zero)
            }
            
            # Save coverage values.
            save(cov_rle, file=paste0(private$filename, ".coverage.RData"))
            
            # Normalize coverages into RPMs.
            rpm_rle = RleList(lapply(cov_rle, "*", 1000000 / private$read_count))
            save(rpm_rle, file=paste0(private$filename, ".coverage_rpm.RData"))
            
            # Save regions.
            cov_regions = reduce(private$all_ranges)
            save(cov_regions, file=paste0(private$filename, ".regions.RData"))
        }
    ),
    private=list(
        filename=NULL,
        con=NULL,
        coverage=data.frame(seqnames=character(0), start=integer(0), end=integer(0), raw=integer(0)),
        all_ranges=GRanges(),
        read_count=0,
        align_line=read_file("input/alignment.txt"),
        chr_size=NULL
    )
)
# Usage examples.
# fake_align_1 = FakeSAM$new("fake_align1.sam")
# fake_align_1$add_flat_coverage("chr1", 1000000, 1005000, 1)
# fake_align_1$write_file_and_close()
# 
# fake_align_1 = FakeSAM$new("fake_align2.sam")
# fake_align_1$add_flat_coverage("chr1", 1000000, 1005000, 2)
# fake_align_1$write_file_and_close()
# 
# fake_align_1 = FakeSAM$new("fake_align3.sam")
# fake_align_1$add_flat_coverage("chr1", 1002500, 1005000, 3)
# fake_align_1$write_file_and_close()

# Example bash commands to go from sam to indexed bam.
# for i in fake_align*.sam
# do
#     bn=`basename $i .sam`
#     samtools view -h -b $i | samtools sort - > $bn.bam
#     samtools index $bn.bam
# done

