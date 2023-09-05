# Linux Commands Cheat Sheet

This commands cheat sheet is particularly for working with and handling **genomic data and files**. Genomic data can be large and convoluted. Open-source bioinformatics programs sometimes require very particular file formats, or sometimes, you simply have 600 sequence files to deal with >.<" Hope this is useful!

## Files

### File Extensions

You will come across MANY different file extensions in bioinformatics. File extensions are arbitary and you can use anything but they should be informative about the data inside them. e.g. I will often name files with a single list inside "something.list" so that I know it's a list just by looking at the file name.

| Ext.            | Name                          | Meaning                                                      |
| --------------- | ----------------------------- | ------------------------------------------------------------ |
| .tar.gz         | zipped tar archive            | Very large data or program inside.                           |
| .fa<br />.fasta | fasta file (generic)          | contains sequence data with a sequence header like: ">seq1". at this stage, a fasta is assembled reads |
| .fna            | fasta file (nucleotide)       | contains nucleotide data (ATGC...)                           |
| .faa            | fasta file (amino acid)       | contains peptide or protein sequences (MSQL...)              |
| .fastq          | unprocessed reads file (nucl) | contains header, sequence and read quality info. Each read is 4 lines of info |
| .tsv            | tab-separated values          | column1    column2    column3                                |
| .csv            | comma-separated values        | column1,column2,column3                                      |
| .gbk            | genbank file                  | has gene/protein specific annotation information. NCBI annotation files |
| .gff            | general feature format        | has gene/protein specific annotation information per contig/chromosome, tab separated |
| sam, bam        | mapping files                 | sam is human-readable mapping file<br />bam is the indexed binary<br />Use [samtools](http://www.htslib.org/doc/samtools.html) or bamtools to deal with these |
| .sh             | shell script                  | has commands to run consecutively using bash                 |

#### Download from internet and unzip Files

Genomic data/databases and bioinformatics programs often come zipped in a tar. Tar = tape archive To [extract](https://linuxize.com/post/how-to-extract-unzip-tar-gz-file/):

```sh
#download database from NCBI ftp site - downloads to current dir
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.faa.tar.gz
file: all.faa.tar.gz

#to extract zipped tar files (.tar.gz)
$ tar -zvxf all.faa.tar.gz 	#-z unzip, -v verbose, -x extract, -f file

#extract to particular location
$ tar -zvxf all.faa.tar.gz -C ./Viral_genes/

#only zipped, not tar
$ gunzip File1.txt.gz
```

### Find files

##### Listing files

```shell
#List the current directory
$ ls
```

```sh
#list a random directory
$ ls ./random/
```

```sh
#list the current directory in chronological order with ownership details
$ ls -lrt
-rw-rw-r-- 1 vmkhot vmkhot  12462377 Feb 16 04:00 viral-fullseq-trim.fa
```

```sh
#list ALL files (including hidden files that start with .)
$ ls -a
.  ..  .bash_history  .bashrc  .cpan  .gnupg  .profile  .ssh  .viminfo
```

```sh
#list in a list
$ ls -1
data
viral-fullseq-trim.fa
gtdbtk_test
```

#### Find and Execute

The "[find](https://man7.org/linux/man-pages/man1/find.1.html)" command can help you search recursively through **directories and subdirectories** for particular file names, extensions or patterns

The general format is:

```sh
find [options] path/to/start expression
```

The simplest type of find command:

```sh
#finds all files with the file extension ".fa" in the random directory and it's subdirectories
$ find ./random -name "*.fa"
```

Other useful find commands/options

```sh
#iname makes the search case-insensitive
$ find ./random -iname "*.fa"

# -regex finding regex patterns
$ find -regex ".*\.fa"

# regex allows you to search multiple extensions; e.g. if you didn't know whether your file was called fasta or fa
$ find -regex ".*\.fa" || ".fasta"

# find FILES with pattern* in random directory, print file name and \n line break
$ find ./random -name "pattern*.fa" -type f printf "%p\n" 

# find and copy all files to new location
find ./ -iname "slurm*.sh" -exec cp {} ./scripts/ \;
```

#### File movement

These are basic file utilities

##### **Copy**

```sh
cp ./location/to/file1 ./new/location/file1
```

##### **Move** or **rename**

```sh
mv ./current/file1 ./destination/file1
```

**Renaming** a file is the same "mv" command, but the current and destination locations are the same directory.

E.g. To rename all ".fa" files to ".fna"

```sh
for f in *.fa; do mv -- "$f" "${f%.fa}.fna"; done
```

To split large directory into subdirectories. Change `300` to wanted number

```bash
i=0; for f in *; do d=dir_$(printf %03d $((i/300+1))); mkdir -p $d; mv "$f" $d; let i++; done
```

##### Remove

!!**<u>STOP</u>**!! Make sure you are absolutely sure about what you are deleting. Linux does not have a "recycle bin". Deleted files and directories are gone FOREVER. 

```sh
#basic
$ rm file.txt

#Delete an empty directory
$ rm -d /path/to/directory

#PARTICULAR DANGER!!!!
#Delete a directory with files and other directories. This will wipe "random" and everything random contains recursively
$ rm -rf path/to/random
```

##### **Move a file to your local desktop**

There's an "scp" command that you can use. IMHO it's annoying to type in all the time. 

I use [MobaXterm](https://mobaxterm.mobatek.net/features.html ) on Windows, which has a terminal, text editor and SFTP connection (Secure File Transfer Protocol), where you can drag and drop files in and out to your desktop. Alternatively you can  [WinSCP](https://winscp.net/eng/docs/introduction#features ) (also Windows, also SFTP). 

For Mac users: [CyberDuck](https://cyberduck.io/sftp/)

### Look inside files

##### **To peek**

```sh
$ less file.txt

$ more file.txt

#To EXIT the less/more view
q

#head will print you FIRST 10 lines from file.txt
$ head file.txt
$ head -n 100 file.txt #prints first 100 lines instead; n = number

#tail will print you LAST 10 files from file.txt
$ tail file.txt
$ tail -n 100 file.txt
```

##### **To edit**

use [nano](https://www.nano-editor.org/) or [vim](https://github.com/vim/vim) - I prefer vim but it's more annoying to use/learn for sure

You can also use WinSCP or MobaXterm  which have in-built text editors.

### Search Inside and Edit Files

Very often you can search inside and edit files right from the command line, without opening the file at all. This is very unintuitive for Windows/Mac users who are used to seeing the changes as they edit, like in Word. 

You might be wondering why one would ever do this. Well, imagine a sequence file with a 1000 sequences and you need to change every header starting with "bin03". This would be painful to do manually, but can be done easily and in seconds with a sed one-liner.  You are likely to encounter these kind of scenarios all the time!

**Grep, awk** and **sed** are the champions of file editing from the command line. They each have lots of options - I've only included here what I've use commonly. 

**Regular Expressions** ([regex](https://tldp.org/LDP/abs/html/x17129.html)) are a way to describe a pattern. Worth learning for quick changes. They are compatible across multiple scripting languages like python and perl as well. Google for more info and how to use these. You can test out your regex patterns here: [Regex101](https://regex101.com/)

#### Grep - Pattern find

You can print lines from a file if you know a pattern you are looking for using "grep" 

[Grep](https://www.gnu.org/software/grep/manual/grep.html#Command_002dline-Options ) is very useful and has tonnes of options (refer to manual) and I use it basically everyday. Grep is pattern matcher so it uses regex patterns. Useful grep one-liners:

```sh
#basic usage
$ grep ">" SeqFile.faa #will print all lines which have > in them (e.g. fasta headers)

#prints all lines with > in them + 2 trailing lines after (fasta header + 2 sequence lines)
$ grep -A 2 ">" SeqFile.faa 
$ grep -B 2    # 2 lines before pattern
$ grep -C 2    # 2 lines before and after pattern

$ grep -v "xxx" #print lines NOT matching xxx

$ grep -i #case insensitive

$ grep -c #print number of lines matching

$ grep -n #print line numbers

$ grep -w #exact match (not regex)

#match multiple patterns
$ grep -A1 "VC37\|VC38\|VC7\|VC36" genes.faa > top4_genes.faa

#My favourite
$ grep -Ff patterns.txt file.txt #this grabs the patterns from a file with a list of patterns and searches for them inside file.txt

#you can combine options
#this will look for headers specified in "headers.list" in an amino acid fasta file "file.faa" and print you the matching headers and 1 sequence line after it
$ grep -A1 -Ff headers.list file.faa
```

To search large **zipped** files without opening them!

```sh
$ zgrep   #zgrep can be used with all above commands
```

#### Awk - Pattern find

[Awk](https://www.gnu.org/software/gawk/manual/gawk.html#Getting-Started ) is effectively a programming language. It's used to search for patterns (line by line) and perform some action on them. It's mega useful for fasta files as well as tsv/csv types.

Basic usage:

This will search for "pattern" in your file and print the entire line "$0" from file.txt

```sh
awk /pattern/ '{print $0}' file.txt
```

##### CSV/TSV 

For the "for loops", see the section on loops.

```shell
# This will split your file by the field separator "-F" (whatever you choose) and print column 1
#tab-separated
$ awk -F'\t' '{print $1}' file.tsv
# comma separated
$ awk -F',' '{print $1}' file.csv

# if statements. Prints full line(row) where column 1 >= 100 in a tab-separated file to output.tsv
$ awk -F'\t' '{if($1 >= 100) print $0{}' file.tsv > output.tsv

# count highest number of columns
$ awk '{print NF}' file | sort -nu | tail -n 1

# AWK can be used to filter tsv files quite effectively. e.g. outputs from blast
# filter column3 (% id) >= 95 AND column11 (evalue) <= 1E-5 and print the full row
$ awk -F'\t' '{if($3 >= 95 && $11 <= 1E-5) print $0}' blast.tsv > filteredblast.tsv

# search strings from filter.txt against column 2 in data.tsv
$ awk -F "\t" 'FNR==NR {hash[$0]; next} !($2 in hash)' filter.txt data.tsv > output.tsv

awk 'BEGIN{FS=OFS="\t"} $1~/VC/ {gsub(/_/, "\|", $1)} 1' temp.out

# filter by column 11 > 0.8, print file name to the first column and cat all files
$ for f in *deduped.tblout; do file=$(basename $f .tblout);awk -v a=${file} '{if ($11 >= 0.8) print a,'\t',$0;}' < ${file}.tblout >> cat_file.tblout ; done
```

##### Fasta

Useful one-liners taken from the internet so not going to explain these

```sh
#SPLIT MULTILINE FASTA INTO INDIVIDUAL FILES
$ awk -F '|' '/^>/ {F=sprintf("%s.fasta",$2); print > F;next;} {print >> F;}' < cyanoWGS_CRISPRs.fasta

#MULTILINE FASTA TO SINGLELINE FASTA
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < input.faa > output.singleline.faa

#SPLIT GIANT FASTA INTO PIECES OF 1000 SEQUENCES (ONLY FOR SINGLELINE)
$ awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < sequences.fa

#UPDATE FASTA HEADERS WITH FILE NAME
$ for f in *.fasta; do file=$(basename $f .fasta); awk -v a=${file} '/^>/{print ">" a "." ++i ; next}{print}' < ${file}.fasta > ${file}.new.fasta; done

```

####  Sed - Search and Replace

"Stream Editor"

I think of [sed](https://www.gnu.org/software/sed/manual/sed.html) primarily as a "search and replace" utility. There are many file changes one can make with sed and they require learning regex and it's very fast for larger files. (much faster than grep imho). 

Commonly used sed one-liners:

```shell
#Basic. the "g" at the end means "global", i.e. replace every instance in the file
sed 's/find/replace/g' input.txt > output.txt

#inplace file editing (!!!DANGER!!!, make sure your command is perfect before using "-i")
sed -i 's/find/replace/g' input.txt

#replace entire first line with xx
sed -i '1s/.*/1620/'

#replace last "_" (underscore) in each line with xx
sed 's/\(.*\)_/\1xx/'

#print lines between two patterns
sed -n '/>pattern1/,/>pattern2/p' file.fasta

# sed if line starts with ">VC"
sed '/^>VC/s/search_string/replace_string/'

# delete line matching pattern (prints the remainder)
sed '/pattern/d'

# delete line matching pattern (exact match with \b) + the line after
# e.g. singleline fasta - it will match the header and delete both header and sequence
sed '/k141_57368\b/,+1 d' file.fa
#OR
sed -e '/k141_57368\b/{N;d;}' file.fa
#multiple
sed '/k141_57368\b\|k141_88513\b/,+1 d' file.fa
```

### Other File Utilities

#### Concatenate

To join files by rows. E.g. concatenate a bunch of fastas into 1 giant fasta

```sh
$ cat file1.fa file2.fa file3.fa > allfiles.fa
```

#### wc -l 

To count the total number of lines

```sh
$ wc -l file.tsv

#multiple files
$ wc -l *.tsv
file1.tsv: 80
file2.tsv: 800
```

Combine with other commands

```sh
#count number of hits with scores above 1000
$ awk '{if ($12 >= 1000) print $0}' blastoutput.tsv | wc -l 
```

#### Sort

Sort will help you sort a list or .tsv or .csv by particular column. Sort is incredibly useful and probably one of my most used commands. I typically use sort to look get a quick idea of my data, e.g. blast outputs.

```sh
#basic and writes to output.txt
$ sort file.txt > output.txt
Apples
Beaches
Cats

#reverse order
$ sort -r file.txt
Cats
Beaches
Apples

#data is numerical
$ sort -n file.txt
10
20
30

#by 2nd column number
$ sort -k 2n file.txt
Apples	100
Cats	200
Beaches	300

#You can also use sort to first sort then deduplicate your rows
$ sort -u file.txt  
```

Sort the top blast hits by bitscore (12th column in outfmt 6).

```sh
#this reverse sorts by 12th column (bitscore) of a blast output file and prints top 10 rows 
$ sort -rk 12n blastoutput.tsv | head 10
```

Or I use it in a pipe from awk to filter data. E.g. If I want to filter all bitscores > 1000 and have them reverse sorted for me

```sh
$ awk '{if ($12 >= 1000) print $0}' blastoutput.tsv | sort -rk 12n > filteredblastout.tsv
```

#### Uniq

Sort is often used with uniq to count or remove duplicates like so: 

```sh
$ sort file.txt | uniq
```

The above command is the same as "sort -u" but uniq has more options too. It's almost always used with the sort command because it sees the "first" instance and doesn't work on unsorted data. 

Uniq is a specialty utility that I don't use very much, I prefer "sort -u".

```sh
#count unique lines
$ sort file.txt | uniq -c

#count duplicates
$ sort file.txt | uniq -d #prints first instance
$ sort file.txt | uniq -D #prints all instances

#remove all duplicates
$ sort file.txt | uniq -u

-i : ignore case
```

#### Cut

[Cut](https://www.geeksforgeeks.org/cut-command-linux-examples/) is useful for extracting data from tab-separated files. Cut requires options to be useful

```sh
$ cut file.tsv	#this produces an error

# extract by bytes
$ cut -b 1,10	#1 - 10 bytes

# extract by character
$ cut -c 1,10	#1 - 10 characters

#extract by field (-f) and delimited (-d)
$ cut -d ' ' -f 2	#delimit file by space and print 2nd field

#cut in action with other commands to get info from headers
input: >contig1_gene1_annotation
$ grep ">" file.fa | cut -d '_' -f 1	
output: >contig1
$ grep ">" file.fa | cut -d '_' -f 1,2,3 --output-delimiter='	'
output: >contig1	gene1	annotation
```

#### Diff

stands for "difference". 

[Diff](https://www.geeksforgeeks.org/diff-command-linux-examples/ ) tells you which lines have to be changed for the files to become identical. Useful for comparing seemingly identical fasta files or blast outputs. 

Diff outputs 3 potential messages.
a: add	c: change	d: delete

```sh
$ diff File1.tsv File2.tsv
```

I don't use this command regularly, I find it easier to compare files by counting rows instead.

### Loops

#### For loop

Sometimes you have to run the same command on many files and you don't want to type it out 100x. You can use a "for loop" to run something on all your files. 

You would put the following code into a shell script file (e.g. do_something.sh) and run in via (bash do_something.sh) Example to adapt:

```sh
for fn in ./*.fa					# iterate through all .fa files in current directory
do							
    echo $fn							# print fasta file name
    newname=$(basename $fn .fa)			# assign name w/o the extension to variable "newname" 
    echo $newname						# print fasta file name
    sample="${newname:1:-6}" 			#"sample" = particular characters of "newname" 
										#(-6 is from right to left)
    echo $sample
done

output:
Sample1_bin11.fa	# $fn
Sample1_bin11		# $newname
Sample1			# $sample
```

Then run by doing:

```sh
$ bash do_something.sh
```

The same loop can also be written as a one-liner

```sh
for fn in ./*.fa; do echo $fn; newname=$(basename $fn .fa); echo $newname; done
```

An actual example to run a mapping command using a for loop

```sh
for fn in ../01_qc/qc_files/CAT*_R1_001.qc.fastq;
do
    base="${fn:0:-16}"
    newname=$(basename $fn .qc.fastq)
    sample="${newname:4:-7}"
    bbmap.sh ref= in=${base}_R1_001.qc.fastq in2=${base}_R2_001.qc.fastq out=./${sample}_bbmap.sam covstats=./${sample}_bbmap_covstats.txt scafstats=./${sample}_bbmap_scafstats.txt threads=20 minid=0.95 ambiguous=toss
done
```
#### While loop

Slightly less useful while loop to loop over the contents of a file

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sh
while IFS= read -r line; do echo "$line"; done < file.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
