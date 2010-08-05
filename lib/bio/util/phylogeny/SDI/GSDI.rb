# = util/phylogeny/SDI/GSDI.rb - generalized SDI 
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an implementation of an algorithm to infer gene duplication
# and speciation events in a rooted gene tree where a node can have more than two children. 
# It is based on the GSDI implementation in forester which is not currently verified. 
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
require 'bio/util/phylogeny/SDI/SDI'

module Bio
  module Algorithm
    # == Description
    #
    # Bio::Algorithm::GSDI implements Speciation/Duplication inferencing for PhyloXML trees. The gene tree can have more than 2 children for
    # a given node.
    #
    # == Usage
    # 
    # Inputs to the constructor are a rooted phyloXML formatted gene tree and a rooted phyloXML formatted binary species tree.
    # The external nodes of the species tree must contain all of the species in the external nodes
    # of the gene tree. 
    #
    # Example:
    # 
    #   require 'bio/util/SDI/GSDI'
    #   # Create an instance of gsdi algorithm
    #   # gene_tree and species_tree are Bio::PhyloXML::Tree objects (see phyloxml documentation)
    #   sdi = Bio::Algorithm::GSDI.new(gene_tree, species_tree)  
    #
    #   # compute the speciation and duplication events:
    #   updated_tree = sdi.compute_speciation_duplications 
    #
    #   # Print number of duplications
    #   # puts sdi.duplications_sum
    #
    class GSDI < SDI
      #hash mapping from species tree node to integer count
      attr_accessor :traversal_counts
      
      # Constructor, calls SDI constructor
      # Create a new Bio::Algorithm::GSDI object
      #   
      #    sdi = Bio::Algorithm::GSDI.new(gene_tree, species_tree)
      #    #gene_tree and species_tree are Bio::PhyloXML::Tree objects (see phyloxml documentation)
      #   
      # The initialization verifies that trees are non-empty and rooted, and that the species tree leaves
      # include all of the species represented in the gene tree and raises an exception if any of the 
      # above conditions fail.
      # ---
        # *Arguments*:
      # * (required) _gene_tree_: Bio::PhyloXML::Tree object
      # * (required) _species_tree_: Bio::PhyloXML::Tree object
      #

      def initialize(gene_tree, species_tree)
        super
        @traversal_counts = {}
      end #initiatilize
      
      # Interprets mapping as determined by LCA method.
      # ---
      # *Arguments*:
      # * (required) _species_node_: Bio::PhyloXML::Tree::Node object
      # * (required) _gene_node_: Bio::PhyloXML::Tree::Node object
      #
      def determine_event_from_mapping(species_node, gene_node)
        sum_g_children_mapping_to_s = 0
        @gene_tree.children(gene_node).each{ |g|
          if @gene_mapping[g] == @species_numbering[species_node]
          sum_g_children_mapping_to_s += 1
          end #if
        }
        # determine sum of traversals
        traversals_sum = 0
        max_traversals = 0
        max_traversals_node = nil
        
        if (!@species_tree.leaves.include?(species_node)) && (species_node != nil)
          spec_node_children = []
        if species_node == @species_tree.root
          spec_node_children = @species_tree.children(@species_tree.root)
        else
          spec_node_children = @species_tree.children(species_node)
        end #if
        spec_node_children.each{ |s|
          if @traversal_counts[s] !=nil
            traversals = @traversal_counts[s]
          else
            traversals = 0
          end
          traversals_sum += traversals
          if traversals > max_traversals
            max_traversals = traversals
          max_traversals_node = s
          end #if
        }
        end #if
        
        if sum_g_children_mapping_to_s > 0
          if traversals_sum == 2
            if gene_node.events == nil
              gene_node.events = Bio::PhyloXML::Events.new
            end #if
            gene_node.events.duplications = 1
            @duplications_sum += 1
          elsif traversals_sum > 2
            if max_traversals <= 1
              if gene_node.events == nil
                gene_node.events = Bio::PhyloXML::Events.new
              end #if
              gene_node.events.speciations = 1
            else
              if gene_node.events == nil
                gene_node.events = Bio::PhyloXML::Events.new
              end #if
              gene_node.events.duplications = 1
              @duplications_sum += 1
              @traversal_counts[max_traversals_node] = 1
            end #if
          else
            if gene_node.events == nil
              gene_node.events = Bio::PhyloXML::Events.new
            end #if
            gene_node.events.duplications = 1
            @duplications_sum += 1
          end #if
        else
          if gene_node.events == nil
            gene_node.events = Bio::PhyloXML::Events.new
          end #if
          gene_node.events.speciations = 1
        end #if
      end #determine_event_from_mapping
      
      # Overrides Bio::Algorithm::SDI::compute_mapping to include non-binary nodes
      # ---
      # *Returns*:: Bio::PhyloXML:Tree object
      #
      def compute_mapping
        root = @gene_tree.root
        if !@gene_tree.leaves.include?(root)
          next_nodes = @gene_tree.children(root)
          next_nodes.each{ |n|
            _compute_mapping(n)
          }
      
          node_mapped_values = []
          next_nodes.each{ |n|
            node_mapped_values.push(@gene_mapping[n])
          }
          max_i = node_mapped_values[0]
          min_i = node_mapped_values[0]
          node_mapped_values.each{ |v|
            if v > max_i
              max_i = v
            end
            if v < min_i
              min_i = v
            end
          }
          while min_i != max_i
            if @traversal_counts[@spec_node_map[max_i]] == nil
              @traversal_counts[@spec_node_map[max_i]] = 1
            else
              @traversal_counts[@spec_node_map[max_i]] += 1
            end #if
            node_mapped_values[node_mapped_values.index(max_i)] = @species_numbering[@species_tree.parent(@spec_node_map[max_i])]
            max_i = node_mapped_values[0]
            min_i = node_mapped_values[0]
            node_mapped_values.each{ |v|
              if v > max_i
                max_i = v
              end
              if v < min_i
                min_i = v
              end
            }  
          end #while
          @gene_mapping[root] = node_mapped_values.last
          determine_event_from_mapping(@spec_node_map[node_mapped_values.last], root)
        end #if
        return @gene_tree
      end #compute_mapping
      
      # Overrides Bio::Algorithm::SDI::_compute_mapping(node)
      # ---
      # *Arguments*:
      # * (required) _node_: Bio::PhyloXML::Node
      #
      def _compute_mapping(node)
      
        if !@gene_tree.leaves.include?(node)
          next_nodes = @gene_tree.children(node)
          next_nodes.each{ |n|
            _compute_mapping(n)
          }
          
          node_mapped_values = []
          next_nodes.each{ |n|
            node_mapped_values.push(@gene_mapping[n])
          }
          max_i = node_mapped_values[0]
          min_i = node_mapped_values[0]
          node_mapped_values.each{ |v|
            if v > max_i
              max_i = v
            end
            if v < min_i
              min_i = v
            end
          }
          while min_i != max_i
            if @traversal_counts[@spec_node_map[max_i]] == nil
              @traversal_counts[@spec_node_map[max_i]] = 1
            else
              @traversal_counts[@spec_node_map[max_i]] += 1
            end #if
            node_mapped_values[node_mapped_values.index(max_i)] = @species_numbering[@species_tree.parent(@spec_node_map[max_i])]
            max_i = node_mapped_values[0]
            min_i = node_mapped_values[0]
            node_mapped_values.each{ |v|
              if v > max_i
                max_i = v
              end
              if v < min_i
                min_i = v
              end
            }
          end #while
          @gene_mapping[node] = node_mapped_values.last
          determine_event_from_mapping(@spec_node_map[node_mapped_values.last], node)
        end #if
      end #_compute_mapping
        
    end #class
  end #module
end #module