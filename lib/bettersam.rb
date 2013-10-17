
class BetterSam

  # meanings of SAM flag components, with index i
  # being one more than the exponent 2 must be raised to to get the
  # value (i.e. value = 2^(i+1))
  $flags = [
    nil,
    0x1,  #    1. read paired 
    0x2,  #    2. read mapped in proper pair (i.e. with acceptable insert size)
    0x4,  #    3. read unmapped
    0x8,  #    4. mate unmapped
    0x10,  #   5. read reverse strand
    0x20,  #   6. mate reverse strand
    0x40,  #   7. first in pair
    0x80,  #   8. second in pair
    0x100,  #  9. not primary alignment
    0x200,  #  10. read fails platform/vendor quality checks
    0x400]  #  11. read is PCR or optical duplicate

public
  attr_accessor :name, :flag, :chrom, :pos, :mapq, :cigar, :mchrom, :mpos, :insert, :seq, :qual, :tags
  attr_accessor :snp
  attr_reader :cigar_list

  def initialize(line=nil)
    # @tags = {}
    parse_line(line) unless line.nil?
  end

  def parse_line(line)
    return false if line[0] == "@"
    
    f = line.chomp.split("\t", -1)
    # raise "SAM lines must have at least 11 fields (had #{f.size})" if f.size < 11

    # colnames = %w(1:name 2:flag 3:chr 4:pos 5:mapq 6:cigar 7:mchr 8:mpos 9:insrt 10:seq 11:qual)

    @name = f[0]
    # @flag = int_or_raw(f[1])
    @flag = f[1].to_i
    @chrom = f[2]
    # @pos = int_or_neg1(f[3])
    @pos = f[3].to_i
    # @mapq = int_or_neg1(f[4])
    # @mapq = f[4]
    @cigar = f[5]
    # @mchrom = f[6]
    # @mpos = int_or_neg1(f[7])
    # @insert = int_or_raw(f[8])
    @seq = f[9]
    @qual = f[10]

    # @tags = {}
    # i = 11
    # while i < f.size
    #   tag = f[i]
    #   i += 1
    #   colon_index = tag.rindex(':') 
    #   raise line if f.rindex == nil
    #   key = tag[0, colon_index]
    #   value = int_or_raw(tag[colon_index + 1, tag.size - colon_index] || "")
    #   @tags[key] = value
    # end

    return true;
  end

  # flag parsing convenience methods

  def read_paired?
    @flag & $flags[1] != 0
  end

  def read_properly_paired?
    @flag & $flags[2] != 0
  end

  def read_unmapped?
    @flag & $flags[3] != 0
  end

  def read_mapped?
    @flag & $flags[3] == 0
  end

  def mate_unmapped?
    @flag & $flags[4] != 0
  end

  def read_reverse_strand?
    @flag & $flags[5] != 0
  end

  def mate_reverse_strand?
    @flag & $flags[6] != 0
  end

  def first_in_pair?
    @flag & $flags[7] != 0
  end

  def second_in_pair?
    @flag & $flags[8] !=0
  end

  def primary_aln?
    !(@flag & $flags[9]) != 0
  end

  def quality_fail?
    @flag & $flags[10] != 0
  end

  def pcr_duplicate?
    @flag & $flags[11] != 0
  end

  # pair convenience methods

  def both_mapped?
    !(self.read_unmapped? && self.mate_unmapped?)
  end

  def pair_opposite_strands?
    (!self.read_reverse_strand? && self.mate_reverse_strand?) || 
      (self.read_reverse_strand? && !self.mate_reverse_strand?)
  end

  def pair_same_strand?
    !self.pair_opposite_strands?
  end

  # cigar parsing methods

  def exact_match?
    @cigar=="100M"
  end

  def endpos
    if !@cigar_list
      self.parse_cigar
    end
    e = @pos
    @cigar_list.each do |h|
      a = h.to_a
      bases = a[0][0]
      match = a[0][1]
      if match =~ /[MD]/
        e += bases
      elsif match == "I"
        e -= bases
      end
    end
    return e
  end

  def parse_cigar
    str = @cigar
    l = str.length
    @cigar_list = []
    while str.length>0
      if str =~ /([0-9]+[MIDNSHPX=]+)/        
        @cigar_list << {$1[0..-2].to_i => $1[-1]}
        str = str.slice($1.length, l)
      else
        puts str
      end
    end
  end

  # snp storing

  def contains_snp?(snp)
    snp >= @pos and snp < self.endpos
  end

  def mark_snp(snp) # snp coordinate is absolute
    if self.contains_snp?(snp)
      if !@cigar_list
        self.parse_cigar
      end
      p = @pos
      s = snp
      @cigar_list.each do |h|
        if p > s and s >= @pos
          @snp = s - @pos
        else
          a = h.to_a
          bases = a[0][0]
          match = a[0][1]
          if match == "M"
            p += bases
          elsif match == "I"
            s += bases
          elsif match == "D"
            s -= bases
          end
        end
      end
      if p > s and s >= @pos
        @snp = s - @pos
      end
    end
    @snp
  end

  def transfer_snp(bs) # load in another bettersam object
    if !self.read_unmapped? and !bs.read_unmapped?
      if bs.seq==nil
        abort "sequence is nil in bs #{bs.name}"
      end
      if bs.snp==nil
        abort "snp is nil in bs #{bs.name}"
      end
      if (self.read_reverse_strand? and bs.read_reverse_strand?) or (!self.read_reverse_strand? and !bs.read_reverse_strand?)
        @snp = bs.snp
      else
        @snp = bs.seq.length - bs.snp
      end
    end
  end

  def put_snp # find the location of a snp on the genome
    if @snp
      if !@cigar_list
        self.parse_cigar
      end
      s = @snp
      p = 0
      @cigar_list.each do |h|
        if p > s
          return s+@pos
        else
          a = h.to_a
          bases = a[0][0]
          match = a[0][1]
          if match=="M"
            p += bases
          elsif match=="D"
            s += bases
          elsif match=="I"
            s -= bases
          end
        end
      end
      if p > s
        return s+@pos
      end
    else
      puts "need to run mark_snp and transfer_snp first"
      return nil
    end
    return -1
  end

  def find_snp(list_of_vcfs)
    if list_of_vcfs!=nil
      prevk=-1
      k = 0
      i = 0
      # list_of_vcfs = @vcfHash[@chrom]
      j = list_of_vcfs.length
      while j - i > 0 and k!=prevk
        prevk = k
        k = (i+j)/2
        # puts "i:#{i} j:#{j} k:#{k}"
        if self.contains_snp?(list_of_vcfs[k].pos)
          i = k
          j = k
          found=true
          # puts "found"
        else
          if @pos > list_of_vcfs[k].pos
            i=k
          else
            j=k
          end
        end
      end
      if found
        return k
      else
        return -1
      end
    else
      return -1
    end
  end

  def get_base_at(p)
    @seq[p]
  end

private

  def int_or_neg1(x)
    Integer(x) rescue -1
  end

  def int_or_raw(x)
    Integer(x) rescue x
  end
end