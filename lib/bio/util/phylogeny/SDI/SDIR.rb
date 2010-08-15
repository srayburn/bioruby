# = util/SDI/SDIR.rb - Binary SDI rerootable trees
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an object that finds the rooting of a gene tree that
# minimizes the number of duplications in the tree. 
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
require 'bio/util/phylogeny/SDI/SDI_rerootable'

module Bio
  module Algorithm
    # == Description
    #
    # Bio::Algorithm::SDIR implements rerooting of gene trees that minimizes the number of duplication
    # events in the tree.
    #
    # == Usage
    #
    #    # Create an instance of SDIR object
    #    # gene_tree and species_tree are binary trees in Bio::PhyloXML::Tree objects
    #    # the gene_tree must be rerootable
    #    sdir = Bio::Algorithm::SDIR.new(gene_tree, species_tree)
    #    
    #    # Find trees with rootings that minimize the number of duplications
    #    trees = sdir.root_and_infer
    #
    class SDIR 
      # integer, minimal number of duplications found
      attr_accessor :_min_dup
      # Bio::PhyloXML::Tree object containing the gene tree
      attr_accessor :gene_tree
      # Bio::PhyloXML::Tree object containing the species tree
      attr_accessor :species_tree
      # integer, number of duplications in final tree
      attr_accessor :duplications_sum
      
  
      # Create an instance of an SDIR object
      #    # gene_tree and species_tree are binary trees in Bio::PhyloXML::Tree objects
      #    # the gene_tree must be rerootable
      #    sdir = Bio::Algorithm::SDIR.new(gene_tree, species_tree)
      # ---
      # *Arguments*:
      # * (required) _gene_tree_: Bio::PhyloXML::Tree object
      # * (required) _species_tree_: Bio::PhyloXML::Tree object
      #
      def initialize(gene_tree, species_tree)
        @gene_tree = gene_tree
        @species_tree = species_tree
        @duplications_sum = 0
      end #initialize
      
      # Finds optimal rootings of tree with respect to minimal duplications.
      #    trees = sdir.root_and_infer
      # --- 
      # *Returns*:: Array of Bio::PhyloXML::Tree objects
      def root_and_infer
        branch_list = []
        rooted_trees = []
        prev_root = nil
        prev_root_c1 = nil
        prev_root_c2 = nil
        
        @_min_dup = 1000000000000
        prev_root_was_dup = false
        # The following is a messy deep copy. Would need to implement clone() for PhyloXML tree object.
        g = Marshal::load(Marshal.dump(gene_tree))
        
        if(g.leaves.length <= 1)
          set_rooted(g, true)
          @_min_dup = 0
          tree_array = [g]
          return tree_array
        end #if
        
        g.nodes.each { |n|
          if !g.leaves.include?(n) && (g.children(n).length != 2)
            raise "Gene tree must be completely binary."
          end #if
        } #end g.nodes.each
        species_tree.nodes.each { |n|
          if !species_tree.leaves.include?(n) && (species_tree.children(n).length != 2)
            raise "Species tree must be completely binary."
          end #if
        } #end species_tree.nodes.each
        set_rooted(g, true)
        #reroot(g, g.leaves[0])
        branches = []
        get_branches_preorder(g, g.root, branches)
        print branches
        
        sdi  = Bio::Algorithm::SDI_rerootable.new(g, species_tree)
        g = sdi.compute_speciation_duplications 
        duplications = sdi.duplications_sum
        used_root_placements = [] + g.children(g.root)
        puts duplications
        puts @_min_dup
        if duplications < @_min_dup
          rooted_trees.clear()
          rooted_trees.push(Marshal::load(Marshal.dump(g)))
          @_min_dup = duplications
          @duplications_sum = @_min_dup
        elsif @duplications == @_min_dup
          rooted_trees.push(Marshal::load(Marshal.dump(g)))
        end #if

        branches.each { |b|
          prev_root = g.root
          puts "previous root:"
          puts prev_root.events
          prev_root_c1 = g.children(prev_root)[0]
          prev_root_c2 = g.children(prev_root)[1]
          prev_root_was_dup = duplication?(prev_root)
          reroot(g,b)
          duplications = sdi.update_mapping(prev_root_was_dup, prev_root_c1, prev_root_c2)
          puts duplications
          puts @_min_dup
          if !used_root_placements.include?(b)
            if duplications < @_min_dup
              rooted_trees.clear()
              rooted_trees.push(Marshal::load(Marshal.dump(g)))
              @_min_dup = duplications
              @duplications_sum = @_min_dup
            elsif @duplications == @_min_dup
              rooted_trees.push(Marshal::load(Marshal.dump(g)))
            end #if
          end #if
          used_root_placements.push(b)
        } #end branches.each
        return rooted_trees
        
      end #root_and_infer
      
      # Set Bio::PhyloXML::Tree.rooted to value
      # ---
      # *Arguments*:
      # * (required) _tree_: Bio::PhyloXML::Tree object
      # * (required) _value_: Boolean
      def set_rooted(tree,value)
        tree.rooted = value
      end #set_rooted
      
      # Reroot the tree. Called by root_and_infer. Modifies tree in place
      # ---
      # *Arguments*:
      # * (required) _tree_: Bio::PhyloXML::Tree object
      # * (required) _new_root_: Bio::PhyloXML::Node object
      def reroot(tree, new_root)
        if new_root == tree.root
          return
        end #if
        if tree.parent(new_root) == tree.root
          return
        end #if
        node1 = tree.children(tree.root)[0]
        node2 = tree.children(tree.root)[1]
        node3 = tree.parent(new_root)
        prev_root = tree.root
        if prev_root == node3
          return
        end #if
        tree.remove_node(tree.root)
        tree.add_edge(node2, node1)
        new_root_node = Bio::PhyloXML::Node.new
        tree.add_node(new_root_node)
        tree.root = new_root_node
        tree.add_edge(tree.root, node3)
        tree.add_edge(tree.root, new_root)
        tree.remove_edge(node3, new_root)
      end #reroot(new_root)
      
      # Get branches of rooted tree in preorder fashion
      # Returns a list of nodes, where the branch is 
      # from the parent of the node to the the node itself.
      # Called by root_and_infer
      # ---
      # *Arguments*:
      # * (required) _tree_: Bio::PhyloXML::Tree object
      # *Returns*:: Array of Bio::PhyloXML::Node objects
      #

      
      def get_branches_preorder(tree, node, b_list)
        b_list.push(node)
        c = tree.children(node)
        if !tree.leaves.include?(node)
          get_branches_preorder(tree, c[0], b_list)
          get_branches_preorder(tree, c[1], b_list)
        end #if
      end #get_branches_preorder(tree, node, b_list)
      
      # Is node a duplication?
      # ---
      # *Arguments*:
      # * (required) _node_: Bio::PhyloXML::Node object
      # *Returns*:: Boolean
      #
      def duplication?(node)
        if node.events != nil
          if node.events.duplications != nil
            return (node.events.duplications > 0)
          else
            return false
          end #if
        else 
          return false
        end #if
      end #duplication?
    end #class
  end #module
end #module