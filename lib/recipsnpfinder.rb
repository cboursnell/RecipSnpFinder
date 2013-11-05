#!/usr/bin/env ruby

require 'rubygems'
require 'trollop'
require 'rasem'

require_relative 'vcf'
require_relative 'bettersam'
require_relative 'snpfinder'

opts = Trollop::options do
  version "v0.0.1a"
  opt :origin, "Origin genome", :required => true, :type => String
  opt :target, "Map fastq files to this genome", :required => true, :type => String
  opt :left, "First set of fastq files", :required => true, :type => String
  opt :right, "Second set of fastq files", :required => true, :type => String
  opt :cores, "Number of cores to use for bowtie", :default => 20, :type => :int
  opt :output, "Output file of location of snps", :required => true, :type => String
  opt :test, "Don't actually run anything"
  opt :deletebam, "Delete bam files"
  opt :verbose, "Be verbose"
end

Trollop::die :origin, "must exist"    if !File.exist?(opts[:origin]) if opts[:origin]
Trollop::die :target, "must exist"    if !File.exist?(opts[:target]) if opts[:target]
Trollop::die :left,   "must exist"    if !File.exist?(opts[:left])   if opts[:left]
Trollop::die :right,  "must exist"    if !File.exist?(opts[:right])  if opts[:right]
Trollop::die :output, "mustn't exist" if File.exist?(opts[:output])  if opts[:output]

snpfinder = SnpFinder.new opts.origin, opts.target, opts.left, opts.right, opts.cores, opts.output, opts.test, opts.verbose

snpfinder.align_origin
snpfinder.align_target

if !snpfinder.vcf_file_exists
  snpfinder.snp_call
else
  puts "vcf file already exists" if opts.verbose
end

if !opts.test
  puts "Loading vcf file" if opts.verbose
  snpfinder.load_vcf_file

  target_index = opts.target.split("/").last.split(".").first
  tmp = opts.left.split("/").last.split(".").first
  sam_target = "#{target_index}-#{tmp}.sam"

  origin_index = opts.origin.split("/").last.split(".").first
  tmp = opts.left.split("/").last.split(".").first
  sam_origin = "#{origin_index}-#{tmp}.sam"

  target_dir = opts.target.split("/")[0..-2].join("/")
  origin_dir = opts.origin.split("/")[0..-2].join("/")

  #debugging
  # sam_target = "nivara-pooled-1000k.sam"
  # sam_origin = "sativa-pooled-10m.sam"
  # sam_target = "barthii-pooled-10m.sam"

  target_fh = File.open("#{target_dir}/#{sam_target}", "r")
  origin_fh = File.open("#{origin_dir}/#{sam_origin}", "r")

  target_line1 = target_fh.readline
  while target_line1 =~ /^@/ # skip comments
    target_line1 = target_fh.readline
  end

  origin_line1 = origin_fh.readline
  while origin_line1 =~ /^@/ # skip comments
    origin_line1 = origin_fh.readline
  end

  target_line2 = target_fh.readline
  origin_line2 = origin_fh.readline
  count=0

  # outputHash = Hash.new
  outputList = Array.new
  puts "Finding sam file length" if opts.verbose
  cmd = "wc -l #{target_dir}/#{sam_target}"
  
  lines = `#{cmd}`
  sam_line_count = lines.to_i
  # sam_line_count = 140178899
  puts "#{cmd} = #{sam_line_count}"


  t1 = BetterSam.new # line 1 from the target sam # could this be speed up by not creating
  t2 = BetterSam.new # line 2 from the target sam # a new bettersam object each time
  o1 = BetterSam.new # line 1 from the origin sam # and reusing ones from previous runs?
  o2 = BetterSam.new # line 2 from the origin sam

  # RubyProf.start

  start_time = Time.now
  progress_time=start_time
  puts "Start time is #{start_time}"
  while target_line1!=nil and origin_line1!=nil and target_line2!=nil and origin_line2!=nil
    count+=2
    if count % 1_000_000==0
      print "#{100*count/sam_line_count}%.. "
      previous_time = progress_time
      progress_time = Time.now
      lines_to_go = sam_line_count - count
      time_taken = progress_time - start_time
      time_to_go = time_taken * lines_to_go / count
      puts "time left: #{time_to_go}\teta: #{progress_time+time_to_go}\t#{progress_time-previous_time}" if opts.verbose
    end
    t1.reinit(target_line1) 
    t2.reinit(target_line2) 
    o1.reinit(origin_line1) 
    o2.reinit(origin_line2) 

    if t1.name == t2.name and o1.name == o2.name and t1.name == o1.name # if the lines are correctly paired
      if t2.first_in_pair?
        tmp = t2
        t2 = t1
        t1 = tmp
      end
      if o2.first_in_pair?
        tmp2 = o1
        o1 = o2
        o2 = tmp2
      end # fix order if bowtie2 screwed it up

      # 1st pair of reads
      if t1.read_mapped? and o1.read_mapped? 
        k = t1.find_snp(snpfinder.list_of_snps(t1.chrom))
        if k>=0

          snp_coord = snpfinder.list_of_snps(t1.chrom)[k].pos
          pos = t1.mark_snp(snp_coord)
          o1.transfer_snp(t1)
          origin_coord = o1.put_snp
          # key_read_name = (t1.name+"/1").to_sym
          outputList << {:target_chrom=>t1.chrom, :origin_chrom=>o1.chrom, :target_coord=>snp_coord, :origin_coord=>origin_coord, :k=>k}
        end
      end

      # 2nd pair of reads
      if t2.read_mapped? and o2.read_mapped? 
        k = t2.find_snp(snpfinder.list_of_snps(t2.chrom))
        if k>=0

          snp_coord = snpfinder.list_of_snps(t2.chrom)[k].pos
          pos = t2.mark_snp(snp_coord)
          o2.transfer_snp(t2)
          origin_coord = o2.put_snp
          # key_read_name = (t2.name+"/2").to_sym
          outputList << {:target_chrom=>t2.chrom, :origin_chrom=>o2.chrom, :target_coord=>snp_coord, :origin_coord=>origin_coord, :k=>k}
        end
      end

    else
      abort "lines are out of sync at #{count}"
    end
    target_line1 = target_fh.readline rescue nil
    origin_line1 = origin_fh.readline rescue nil
    target_line2 = target_fh.readline rescue nil
    origin_line2 = origin_fh.readline rescue nil
  end


  # build output
  sum_hash = Hash.new
  # print "Building output..."
  output = ""
  count=0
  start_time2 = Time.now
  progress_time=start_time2
  # size = outputHash.keys.size
  size = outputList.size

  puts "Opening file #{opts.output} for appending" if opts.verbose

  # outputHash.each_pair do |key, value|
  outputList.each do |value|
    count+=1
    # value is a list of hashes
    if count % 1_000_000==0
      previous_time = progress_time
      progress_time = Time.now
      lines_to_go = size - count
      time_taken = progress_time - start_time2
      time_to_go = time_taken * lines_to_go / count
      puts "#{100*count/size}%..\ttime left: #{time_to_go}\teta: #{progress_time+time_to_go}\t#{progress_time-previous_time}" if opts.verbose

    end
    # output << "#{key}\t#{value[:origin_chrom]}\t#{value[:origin_coord]}\t#{value[:target_chrom]}\t#{value[:target_coord]}\t#{value[:k]}\n"
    output << "#{value[:origin_chrom]}\t#{value[:origin_coord]}\t#{value[:target_chrom]}\t#{value[:target_coord]}\t#{value[:k]}\n"
  end

  puts "Done"

  print "Writing output to file..."
  File.open("#{opts.output}", "w") {|file| file.write(output)}
  puts "Done"

  finish_time = Time.now
  puts "Finish time was #{finish_time}"
  puts "Took #{(finish_time-start_time).round(1)}s"
end
