#!/usr/bin/env ruby

require 'rubygems'
require 'trollop'
require 'rasem'
require 'ruby-prof'

require_relative 'vcf'
require_relative 'bettersam'
require_relative 'snpfinder'

#
# i think this takes too long to run for what it is doing. i think i need to rewrite this in java
#

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
snpfinder.create_origin_index
snpfinder.align_origin
snpfinder.create_target_index
snpfinder.align_target
snpfinder.snp_call

if !opts.test
  puts "Loading vcf file" if opts.verbose
  snpfinder.load_vcf_file

  target_index = opts.target.split("/").last.split(".").first
  tmp = opts.left.split("/").last.split(".").first
  sam_target = "#{target_index}-#{tmp}.sam"

  origin_index = opts.origin.split("/").last.split(".").first
  tmp = opts.left.split("/").last.split(".").first
  sam_origin = "#{origin_index}-#{tmp}.sam"

  dir = opts.target.split("/")[0..-2].join("/")

  #debugging
  # sam_target = "nivara-pooled-100k.sam"
  sam_target = "barthii-pooled-10m.sam"
  sam_origin = "sativa-pooled-10m.sam"

  target_fh = File.open("#{dir}/#{sam_target}", "r")
  origin_fh = File.open("#{dir}/#{sam_origin}", "r")

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

  outputHash = Hash.new
  puts "Finding sam file length" if opts.verbose
  cmd = "wc -l #{dir}/#{sam_target}"
  
  lines = `#{cmd}`
  sam_line_count = lines.to_i
  # sam_line_count = 100000
  puts "#{cmd} #{sam_line_count}"

  start_time = Time.now
  puts "Start time is #{start_time}"
  while target_line1!=nil and origin_line1!=nil and target_line2!=nil and origin_line2!=nil
    count+=2
    if count % 100000==0
      print "#{100*count/sam_line_count}%.. "
      # progress = sam_line_count.to_f/count.to_f
      # progress_time = Time.now
      # time_taken = progress_time - start_time
      # eta = progress * time_taken
      # puts "eta: #{start_time+eta}"
      progress_time = Time.now
      lines_to_go = sam_line_count - count
      time_taken = progress_time - start_time
      time_to_go = time_taken * lines_to_go / count
      puts "time left: #{time_to_go}\teta: #{progress_time+time_to_go}"
    end
    t1 = BetterSam.new(target_line1) # line 1 from the target sam
    t2 = BetterSam.new(target_line2) # line 2 from the target sam
    o1 = BetterSam.new(origin_line1) # line 1 from the origin sam
    o2 = BetterSam.new(origin_line2) # line 2 from the origin sam

    if t1.name == t2.name and o1.name == o2.name # if the lines are correctly paired
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
          if !outputHash.has_key?(k)
            outputHash[k]=[]
          end
          snp_coord = snpfinder.list_of_snps(t1.chrom)[k].pos
          pos = t1.mark_snp(snp_coord)
          o1.transfer_snp(t1)
          origin_coord = o1.put_snp
          outputHash[k] << {:name=> t1.name, :target_chrom=>t1.chrom, :origin_chrom=>o1.chrom, :target_coord=>snp_coord, :pos=>pos, :origin_coord=>origin_coord}
        end
      end
      # 2nd pair of reads
      if t2.read_mapped? and o2.read_mapped? 
        k = t2.find_snp(snpfinder.list_of_snps(t2.chrom))
        if k>=0
          if !outputHash.has_key?(k)
            outputHash[k]=[]
          end
          snp_coord = snpfinder.list_of_snps(t2.chrom)[k].pos
          pos = t2.mark_snp(snp_coord)
          o2.transfer_snp(t2)
          origin_coord = o2.put_snp
          outputHash[k] << {:name=> t2.name, :target_chrom=>t2.chrom, :origin_chrom=>o2.chrom, :target_coord=>snp_coord, :pos=>pos, :origin_coord=>origin_coord}
        end
      end
    else
      abort "lines are out of sync"
    end
    target_line1 = target_fh.readline rescue nil
    origin_line1 = origin_fh.readline rescue nil
    target_line2 = target_fh.readline rescue nil
    origin_line2 = origin_fh.readline rescue nil
  end


  RubyProf.start

  # build output
  sum_hash = Hash.new
  print "Building output..."
  output = ""
  count=0
  start_time2 = Time.now
  size = outputHash.keys.size
  outputHash.each_pair do |key, value|
    count+=1
    # value is a list of hashes
    if count % 100==0
      previous_time = progress_time
      progress_time = Time.now
      lines_to_go = size - count
      time_taken = progress_time - start_time2
      time_to_go = time_taken * lines_to_go / count
      puts "#{100*count/size}%..\ttime left: #{time_to_go}\teta: #{progress_time+time_to_go}\t#{progress_time-previous_time}"
    end
    value.each do |snp_hash|
      output += "#{key}\t#{snp_hash[:name]}\t#{snp_hash[:target_chrom]}\t#{snp_hash[:target_coord]}\t#{snp_hash[:origin_chrom]}\t#{snp_hash[:origin_coord]}\t#{snp_hash[:pos]}\n"
      # if !sum_hash.has_key?(snp_hash[:target_chrom])
      #   sum_hash[snp_hash[:target_chrom]] = Hash.new
      # end
      # if !sum_hash[snp_hash[:target_chrom]].has_key?(snp_hash[:target_coord])
      #   sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]] = Hash.new
      # end
      # if !sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]].has_key?(snp_hash[:origin_chrom])
      #   sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]][snp_hash[:origin_chrom]]=Hash.new
      # end
      # if !sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]][snp_hash[:origin_chrom]].has_key?(snp_hash[:origin_coord])
      #   sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]][snp_hash[:origin_chrom]][snp_hash[:origin_coord]]=0
      # end
      # sum_hash[snp_hash[:target_chrom]][snp_hash[:target_coord]][snp_hash[:origin_chrom]][snp_hash[:origin_coord]]+=1
    end
  end

  # sum_output=""
  # sum_hash.each_pair do |key,value|
  #   # puts "target chromosome: #{key}"
  #   value.each_pair do |key2, value2|
  #     # puts "  target coord: #{key2}"
  #     value2.each_pair do |key3, value3|
  #       # puts "    origin chromosome: #{key3}"
  #       value3.each_pair do |key4, value4|
  #         sum_output += "#{key}\t#{key2}\t#{key3}\t#{key4}\t#{value4}\n"
  #       end
  #     end
  #   end
  # end
  puts "Done"

  # write output to file

  print "Writing output to file..."
  File.open("#{opts.output}", "w") {|file| file.write(output)}
  # File.open("sum_#{opts.output}", "w") {|file| file.write(sum_output)}
  puts "Done"

  finish_time = Time.now
  puts "Finish time was #{finish_time}"
  puts "Took #{(finish_time-start_time).round(1)}s"


  result = RubyProf.stop

  # Print a flat profile to text
  printer = RubyProf::FlatPrinter.new(result)
  printer.print(STDOUT)
end
