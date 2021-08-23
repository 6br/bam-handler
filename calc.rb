class Array
    def sum
      reduce(:+)
    end
  
    def mean
      sum.to_f / size
    end
  
    def median
      sorted = self.sort
      len = sorted.length
      return 0 if len == 0
      (sorted[(len - 1) / 2] + sorted[len / 2]) / 2.0
    end
  
    def var
      m = mean
      reduce(0) { |a,b| a + (b - m) ** 2 } / (size - 1)
    end
  
    def sd
      Math.sqrt(var)
    end
  end
  
  bin_range = [1, 9, 73, 585, 4681]
  list_chunk = [[],[],[],[],[]]
  list_offset = [[],[],[],[],[]]
  
  while line = gets
    k = line.chomp.split("\t")
    idx = bin_range.index{|t| t > k[0].to_i}
    idx = 4 unless idx
    list_chunk[idx] << k[1].to_i
    list_offset[idx] << k[2].to_i
  end
  

  puts ["layer", "length", "sum", "mean", "var", "median"].join("\t")
  list_chunk.each.with_index{|t, i| puts [i, t.length, t.sum, t.mean, t.var, t.median].join("\t")}

#  puts ["layer", "length", "sum", "mean", "var", "median"].join("\t")
  list_offset.each.with_index{|t, i| puts [i, t.length, t.sum, t.mean, t.var, t.median].join("\t")}