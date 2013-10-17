class Vcf
  attr_accessor :contig, :pos, :type

  def initialize(readline)
    cols = readline.split("\t")
    @contig = cols[0]
    @pos = cols[1].to_i
    if cols[7]=~ /INDEL/
      @type = "INDEL"
    else
      @type = "SNP"
    end
  end

  def to_s
    "#{@contig} #{@pos} #{@type}"
  end
end