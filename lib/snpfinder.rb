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
      @origin_dir = @origin.split("/")[0..-2].join("/")
    else
      @origin_dir = "."
    end
    if @target.split("/").length > 1
      @target_dir = @target.split("/")[0..-2].join("/")
    else
      @target_dir = "."
    end
  end

  def align_origin
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @origin_index = @origin.split("/").last.split(".").first
    @origin_sam = "#{@origin_index}-#{tmp}.sam"

    if !File.exists?("#{@origin_dir}/#{@origin_sam}")
      if !File.exists?("#{@origin_dir}/#{@origin_index}.1.bt2")
        cmd = "bowtie2-build #{@origin} #{@origin_dir}/#{@origin_index}"
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      else
        puts "#{@origin_dir}/#{@origin_index}.1.bt2 already exists" if @verbose
      end


      cmd = "bowtie2 -t -p #{@cores} -x #{@origin_dir}/#{@origin_index} -1 #{@left} -2 #{@right} --very-sensitive -S #{@origin_dir}/#{@origin_sam} --reorder" 
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@origin_dir}/#{@origin_sam} already exists" if @verbose
    end
  end

  def align_target
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @target_index = @target.split("/").last.split(".").first
    @target_sam = "#{@target_index}-#{tmp}.sam"

    if !File.exists?("#{@target_dir}/#{@target_sam}")
      if !File.exists?("#{@target_dir}/#{@target_index}.1.bt2")
        cmd = "bowtie2-build #{@target} #{@target_dir}/#{@target_index}"
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      else
        puts "#{@target_dir}/#{@target_index}.1.bt2 already exists" if @verbose
      end

      cmd = "bowtie2 -t -p #{@cores} -x #{@target_dir}/#{@target_index} -1 #{@left} -2 #{@right} --very-sensitive -S #{@target_dir}/#{@target_sam} --reorder" 
      puts " $ #{cmd}"
      `#{cmd}` if !@test
    else
      puts "#{@target_dir}/#{@target_sam} already exists" if @verbose
    end
  end

  def snp_call
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    if !File.exists?("#{@target_dir}/#{@target_index}-#{tmp}.vcf")

      bam = "#{@target_sam.split(".").first}.bam"
      if !File.exists?("#{@target_dir}/#{bam}")
        ### make bam file
        cmd = "samtools view -bS -o #{@target_dir}/#{bam} #{@target_dir}/#{@target_sam}"
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      else
        puts "#{bam} file already exists" if @verbose
      end


      ### sort bam file
      sorted = "#{bam.split(".").first}.sorted"
      if !File.exists?("#{@target_dir}/#{sorted}.bam")
        cmd = "samtools sort #{@target_dir}/#{bam} #{@target_dir}/#{sorted}"
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      else
        puts "#{sorted}.bam file already exists" if @verbose
      end

      samtools=false
      if samtools 
        ### do snp calling
        cmd = "samtools mpileup -uf #{@target} #{@target_dir}/#{sorted}.bam | bcftools view -bvcg - > #{@target_dir}/#{@target_index}-#{tmp}.bcf"; 
        puts " $ #{cmd}"
        `#{cmd}` if !@test

        ### convert bcf to vcf
        cmd = "bcftools view #{@target_dir}/#{@target_index}-#{tmp}.bcf | vcfutils.pl varFilter -D100 > #{@target_dir}/#{@target_index}-#{tmp}.vcf";
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      else

        if !File.exists?("#{@target}.fai")
          cmd = "samtools faidx #{@target}"
          puts cmd
          `#{cmd}` if !@test
        else
          puts "#{@target}.fai already exists"
        end

        if !File.exists?("#{@target_dir}/#{sorted}.bam.bai")
          cmd = "samtools index #{@target_dir}/#{sorted}.bam"
          puts cmd
          `#{cmd}` if !@test
        else
          puts "#{@target_dir}/#{sorted}.bam.bai already exists"
        end

        ### platypus
        cmd = "python /home/cmb211/platypus/Platypus.py callVariants --bamFiles=#{@target_dir}/#{sorted}.bam  --output=#{@target_dir}/#{@target_index}-#{tmp}.vcf --refFile=#{@target}"
        puts " $ #{cmd}"
        `#{cmd}` if !@test
      end
    else
      puts "#{@target_dir}/#{@target_index}-#{tmp}.vcf vcf file already exists" if @verbose
    end
  end

  def vcf_file_exists
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @target_index = @target.split("/").last.split(".").first
    @origin_index = @origin.split("/").last.split(".").first
    if File.exists?("#{@target_dir}/#{@target_index}-#{tmp}.vcf")
      return true
    else
      return false
    end
  end

  def load_vcf_file
    tmp = @left.split("/").last.split(".").first # tmp is a name like 'pooled'
    @vcfHash = Hash.new
    File.open("#{@target_dir}/#{@target_index}-#{tmp}.vcf", "r").each_line do |line|
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