# = util/SDI/sdi_binary.rb - Binary SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#
# == Description
#
# This file contains an implementation of an algorithm to infer gene duplication
# and speciation events in a binary gene tree. 
#
# == Usage
# 
# Inputs are a phyloXML formatted binary gene tree and a phyloXML formatted binary species tree.
# The external nodes of the species tree must contain all of the species in the external nodes
# of the gene tree. The output is an updated phyloXML formatted gene tree.
#
# To create an instance of the algorithm:
#
# sdi = SDI.new('species.xml', 'gene.xml')
#
# where 'species.xml' is replaced by the path to the species tree
# and 'gene.xml' is replaced by the path to the gene tree
#
# To compute the speciation and duplication events:
#
# sdi.computeSpeciationDuplications 
#
# and to write the update gene tree output:
#
# sdi.writeUpdatedGeneTreePhyloXML('updated_output.xml')
#
# where 'updated_output.xml' is replaced by the path to the desired output file.
#
# == References
#
# * http://wiki.github.com/srayburn/bioruby
# * http://bioinformatics.oxfordjournals.org/cgi/reprint/17/9/821
# * http://www.phyloxml.org
# * https://www.nescent.org/wg_phyloinformatics/PhyloSoC:PhyloXML_support_in_BioRuby
# 

#  load libraries required

require 'bio/db/phyloxml/phyloxml_parser'

module Bio
class SDI

  # rooted binary gene tree
  attr_accessor :gene_tree
  # rooted binary species tree of all species in gene tree
  attr_accessor :species_tree
  # mapping from node in gene tree to integer number
  attr_accessor :gene_mapping 
  # mapping from node in species tree to integer number
  attr_accessor :species_numbering
  # mapping from integer number to node in species tree
  attr_accessor :spec_node_map

  # The parameters are phyloxml tree objects containing
  # the gene tree and the species tree.
  # Raises exception when either tree is empty, unrooted, or 
  # when the species tree does not contain all species in the gene tree.
  def initialize(gene_tree, species_tree)
    # phyloxml trees 
    @gene_tree = gene_tree
    @species_tree = species_tree

    # Raise exception for empty tree
    if (@gene_tree == nil) || (@species_tree == nil)
      raise "Gene and Species trees must be non-empty."
    end #if

    # raise exception for unrooted tree
    if !isRooted?(@gene_tree) || !isRooted?(@species_tree)
      raise "Gene and Species trees must be rooted."
    end # if
    
    # raise exception if species tree does not contain all of gene tree's species
    if !isSubset?(@species_tree, @gene_tree)
      raise "All species in gene tree leaf level must be represented in species tree."
    end #if
    
    @gene_mapping = {}
    @species_numbering = {}
    @spec_node_map = {}
   
  end #initialize
  
  # Wrapper method for the algorithm. Calls initialization and mapping computation.
  def computeSpeciationDuplications

    initializeMapping
    computeMapping

  end #computeSpeciationDuplications

  # Writes updated gene tree to phyloxml file. filename contains path to new file.
  def writeUpdatedGeneTreePhyloXML(filename)
   
    # Create new phyloxml writer
    writer = Bio::PhyloXML::Writer.new(filename)
   
    # Write tree to the file tree.xml
    writer.write(@gene_tree)

  end #writeUpdatedGeneTreePhyloXML

  # Computes the mapping between gene and species trees and updates each clade in 
  # gene tree to have either a speciation or duplication event. Calls recursive method
  # _computeMapping(node). Called by computeSpeciationDuplications()
  
  def computeMapping 
 
    root = @gene_tree.root
    nextNodes = @gene_tree.children(root)
    nextNodes.each { |n|
      index = _computeMapping(n) 
    }
    # this is how forester does it…
    a = @gene_mapping[nextNodes[0]]
    b = @gene_mapping[nextNodes[1]]
    
    @gene_mapping[root] = @species_numbering[@species_tree.lowest_common_ancestor(@spec_node_map[a], @spec_node_map[b])]

    if (@gene_mapping[root] == @gene_mapping[nextNodes[0]]) || (@gene_mapping[root] == @gene_mapping[nextNodes[1]])
      root.events = Bio::PhyloXML::Events.new
      root.events.duplications = 1
    else
      root.events = Bio::PhyloXML::Events.new
      root.events.speciations = 1
    end #if

  end #computeMapping
  
  
  # Recursive, private helper to computeMapping(). 
  def _computeMapping(node)  		 #:doc:
     
    nextNodes = @gene_tree.children(node)
    nextNodes.each { |n|
      index = _computeMapping(n) 
    }
    if nextNodes[0] != nil:
      a = @gene_mapping[nextNodes[0]]
      if nextNodes[1] != nil:
        b = @gene_mapping[nextNodes[1]]
      else 
        b = @gene_mapping[nextNodes[0]]
      end #if  
    
      @gene_mapping[node] = @species_numbering[@species_tree.lowest_common_ancestor(@spec_node_map[a], @spec_node_map[b])]
      if (@gene_mapping[node] == @gene_mapping[nextNodes[0]]) || (@gene_mapping[node] == @gene_mapping[nextNodes[1]])
        node.events = Bio::PhyloXML::Events.new
        node.events.duplications = 1
      else
        node.events = Bio::PhyloXML::Events.new
        node.events.speciations = 1
      end #if

    end #if

  end #_computeMapping
  
  # Numbers nodes of species tree in a preorder fashion.
  # Calls private helper _initializeSpeciesMap()
  # modifies instance variables @species_numbering and @spec_node_map
  
  def initializeSpeciesMap
    root = @species_tree.root
    index = 1
    
    @species_numbering[root] = index
    @spec_node_map[index] = root
    index += 1
    
    nextNodes = species_tree.children(root)
    nextNodes.each { |n|
      index = _initializeSpeciesMap( n, index)
    }
   
  end #initializeSpeciesMap
 
  # recursive, private helper of initializeSpeciesMap()
  def _initializeSpeciesMap( node, index)         # :doc:
   
    @species_numbering[node] = index
    @spec_node_map[index] = node
    index = index + 1
   
    nextNodes = species_tree.children(node)
    nextNodes.each{ |n|
      index = _initializeSpeciesMap( n, index)
    }
    return index

  end #_initializeSpeciesMap
  
  # Implements initialization of algorithm.
  # Numbers species tree nodes in preorder traversal
  # and sets mapping of leaf nodes in gene tree to
  # the mapping of a matching leaf node in the species tree.
  def initializeMapping
    
    initializeSpeciesMap
    index = 0
    @gene_tree.leaves.each{ |n|
      @species_tree.leaves.each{ |s|
        if nodeEqual?(n, s)
          @gene_mapping[n] = @species_numbering[s]
        end # if
       }
    }
 
  end #initializeMapping

  # Tests to see if nodes are equivalent. Checks taxonomy id first, then code, then scientific name, then common name.
  # Raises fatal exception if nodes do not have enough information to match.
  def nodeEqual?(node1, node2)
    if (!node1.taxonomies.empty?) && (!node2.taxonomies.empty?)
    # compare taxonomy id if exists and provider is the same
    if (node1.taxonomies[0].taxonomy_id != nil) && (node2.taxonomies[0].taxonomy_id != nil)
      if (node1.taxonomies[0].taxonomy_id.provider != nil) && (node2.taxonomies[0].taxonomy_id.provider != nil)
        if node1.taxonomies[0].taxonomy_id.value == node2.taxonomies[0].taxonomy_id.value
         return true
        else 
          return false
        end # if
      end #if
    
    # otherwise compare code
    elsif (node1.taxonomies[0].code != nil) && (node2.taxonomies[0].code != nil) 
      if node1.taxonomies[0].code == node2.taxonomies[0].code
        return true
      else 
        return false
      end #if
  
    # otherwise compare scientific name
    elsif (node1.taxonomies[0].scientific_name != nil) && (node2.taxonomies[0].scientific_name != nil)
      if node1.taxonomies[0].scientific_name == node2.taxonomies[0].scientific_name
        return true
      else 
        return false
      end #if
  
    # otherwise compare common names
    elsif (node1.taxonomies[0].common_names[0] != nil) && (node2.taxonomies[0].common_names[0] != nil)
      if node1.taxonomies[0].common_names[0] == node2.taxonomies[0].common_names[0]
        return true
      else 
        return false
      end # if
    end #if  
    end #if
    # otherwise, not enough in common to compare
    raise "Nodes must share an identifier to be compared."
    

  end #nodeEqual?

  # Tests tree for root, boolean returning.
  def isRooted?(tree)
    if tree.root == nil:
      return false
    else
      return true
    end #if
    
  end #isRooted?
 
  # Tests that all leaves in subset are included in leaves of superset
  def isSubset?(superset, subset)

    superset_leaves = superset.leaves
    subset_leaves = subset.leaves
    
    subset_leaves.each { |node|
      inner = false
      superset_leaves.each{ |node1|
        inner = false
        if nodeEqual?(node, node1)
          inner = true
          break
        end #if
      } #superset_leaves.each
      if !inner
        return false
      end #if
    } #subset_leave.each
    return true

  end #isSubset?
 
  private :_initializeSpeciesMap
  private :_computeMapping

end
end
