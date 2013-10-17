class SnpFinder
  public
  attr_accessor :origin, :target, :left, :right, :output, :test, :verbose

  def initialize(origin, target, left, right, cores, output, test, verbose)
    @origin = origin
    @target = target
    @left = left
    @right = right
    @cores = cores
    @output = output
    @test = test
    @verbose = verbose
    if @origin.split("/").length > 1
      @dir = @origin.split("/")[0..-2].join("/")
    else
      @dir = "."
    end
  end

  def create_origin_index
    @origin_index = @origin.split("/").last.split(".").first
    if !File.exists?("#{@dir}/#{@origin_index}.1.bt2")
      cmd = "bowtie2-build #{@origin} #{@dir}/#{@origin_index}"
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@dir}/#{@origin_index}.1.bt2 already exists" if @verbose
    end
  end

  def create_target_index
    @target_index = @target.split("/").last.split(".").first
    if !File.exists?("#{@dir}/#{@target_index}.1.bt2")
      cmd = "bowtie2-build #{@target} #{@dir}/#{@target_index}"
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@dir}/#{@target_index}.1.bt2 already exists" if @verbose
    end
  end

  def align_origin
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @origin_sam = "#{@origin_index}-#{tmp}.sam"

    if !File.exists?("#{@dir}/#{@origin_sam}")
      cmd = "bowtie2 -t -p #{@cores} -x #{@dir}/#{@origin_index} -1 #{@left} -2 #{@right} --very-sensitive -S #{@dir}/#{@origin_sam} --reorder" 
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@dir}/#{@origin_sam} already exists" if @verbose
    end
  end

  def align_target
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @target_sam = "#{@target_index}-#{tmp}.sam"

    if !File.exists?("#{@dir}/#{@target_sam}")
      cmd = "bowtie2 -t -p #{@cores} -x #{@dir}/#{@target_index} -1 #{@left} -2 #{@right} --very-sensitive -S #{@dir}/#{@target_sam} --reorder" 
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@dir}/#{@target_sam} already exists" if @verbose
    end
  end

  def snp_call
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    if !File.exists?("#{@dir}/#{@target_index}-#{tmp}.vcf")

      ### make bam file
      bam = "#{@target_sam.split(".").first}.bam"
      cmd = "samtools view -bS -o #{@dir}/#{bam} #{@dir}/#{@target_sam}"
      puts " $ #{cmd}"
      `#{cmd}` if !@test

      ### sort bam file
      sorted = "#{bam.split(".").first}.sorted"
      cmd = "samtools sort #{@dir}/#{bam} #{@dir}/#{sorted}"
      puts " $ #{cmd}"
      `#{cmd}` if !@test

      ### do snp calling
      cmd = "samtools mpileup -uf #{@origin} #{@dir}/#{sorted}.bam | bcftools view -bvcg - > #{@dir}/#{@target_index}-#{tmp}.bcf"; 
      puts " $ #{cmd}"
      `#{cmd}` if !@test

      ### convert bcf to vcf
      cmd = "bcftools view #{@dir}/#{@target_index}-#{tmp}.bcf | vcfutils.pl varFilter -D100 > #{@dir}/#{@target_index}-#{tmp}.vcf";
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@dir}/#{@target_index}-#{tmp}.vcf vcf file already exists" if @verbose
    end
  end

  def vcf_file_exists
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @target_index = @target.split("/").last.split(".").first
    @origin_index = @origin.split("/").last.split(".").first
    if File.exists?("#{@dir}/#{@target_index}-#{tmp}.vcf")
      return true
    else
      return false
    end
  end

  def load_vcf_file
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @vcfHash = Hash.new
    File.open("#{@dir}/#{@target_index}-#{tmp}.vcf", "r").each_line do |line|
      if line !~ /^#/
        cols = line.split("\t")
        if !@vcfHash.has_key?(cols[0]) # cols[0] is the conting/chromosome the snp was found in
          @vcfHash[cols[0]] = []
        end
        if cols[3].length == cols[4].length # only store snps not indels; TOO CONFUSING!
          @vcfHash[cols[0]] << Vcf.new(line)
        end
      end
    end
  end

  def list_of_snps(chrom)
    # puts "list_of_snps(#{chrom})"
    if @vcfHash.has_key?(chrom)
      return @vcfHash[chrom]
    else
      # puts "returning nil because I couldn't find #{chrom}"
      return nil
    end
  end
end