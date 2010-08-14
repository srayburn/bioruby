# = util/phylogeny/SDI/SDI_rerootable.rb - Binary SDI rerootable trees
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an implementation of an algorithm to infer gene duplication
# and speciation events in a rooted binary gene tree. 
#
#
#
# == References
#
# * http://wiki.github.com/srayburn/bioruby
# * http://bioinformatics.oxfordjournals.org/cgi/reprint/17/9/821
# * http://www.phyloxml.org
# * https://www.nescent.org/wg_phyloinformatics/PhyloSoC:PhyloXML_support_in_BioRuby
# 

#  load libraries required
#require 'lib/bio/util/phylogeny/SDI/SDI'

module Bio
  module Algorithm
    class SDI_rerootable < SDI
      def update_mapping(prev_root_was_dup, prev_root_c1, prev_root_c2)
        if (@gene_tree.children(@gene_tree.root)[0] == prev_root_c1) || (@gene_tree.children(@gene_tree.root)[1] == prev_root_c1)
          calculate_mapping_for_node( prev_root_c1)
        else
          calculate_mapping_for_node(prev_root_c2)
        end #if
        
        if prev_root_was_dup
          if @gene_tree.root.events == nil
          @gene_tree.root.events = Bio::PhyloXML::Events.new
        end #if
          @gene_tree.root.events.duplications = 1
        else
          if @gene_tree.root.events == nil
          @gene_tree.root.events = Bio::PhyloXML::Events.new
        end #if
          @gene_tree.root.events.speciations = 1
        end #if
        
        calculate_mapping_for_node(@gene_tree.root)
        return @duplications_sum
      end #update_mapping
    
      def calculate_mapping_for_node(node)
        if !@gene_tree.leaves.include?(node)
          if node.events != nil
          if node.events.duplications != nil && node.events.duplications > 0
            was_duplication = true
          end #if
        else 
          was_duplication = false
        end #if
          a = @gene_mapping[@gene_tree.children(node)[0]]
            b = @gene_mapping[@gene_tree.children(node)[1]]
        
            @gene_mapping[node] = @species_numbering[@species_tree.lowest_common_ancestor(@spec_node_map[a], @spec_node_map[b])]
     
            if (@gene_mapping[node] == @gene_mapping[@gene_tree.children(node)[0]]) || (@gene_mapping[node] == @gene_mapping[@gene_tree.children(node)[1]])
              if node.events == nil
          node.events = Bio::PhyloXML::Events.new
          end #if
              node.events.duplications = 1
          if !was_duplication
            @duplications_sum += 1
          end #if
            else
          if node.events == nil
                node.events = Bio::PhyloXML::Events.new
              end #if
          node.events.speciations = 1
            end #if
          end #if
      end  #calculate_mapping_for_node(node)
    end #class
  end #module
end #module