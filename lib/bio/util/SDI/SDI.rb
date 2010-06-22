# = util/SDI/sdi_binary.rb - Binary SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#

#  load libraries required

require 'bio/db/phyloxml/phyloxml_parser'


class SDI

  attr_accessor :gene_tree
  attr_accessor :species_tree
  attr_accessor :gene_mapping 
  attr_accessor :species_numbering
  attr_accessor :spec_node_map
  attr_accessor :gene_node_map

  def initialize(gene_filename, species_filename)
    # load phyloxml trees (assumption: each file contains only 1 tree)
    @gene_tree = Bio::PhyloXML::Parser.open(gene_filename)
    @species_tree = Bio::PhyloXML::Parser.open(species_filename)
    @gene_tree = @gene_tree.next_tree
    @species_tree = @species_tree.next_tree

    # test trees non-empty
    if (@gene_tree == nil) || (@species_tree == nil)
      raise "Gene and Species trees must be non-empty."
    end #if

    # test rooted
    if !isRooted?(@gene_tree) || !isRooted?(@species_tree)
      raise "Gene and Species trees must be rooted."
    end # if
    
    # test subset
    if !isSubset?(@species_tree, @gene_tree)
      raise "All species in gene tree leaf level must be represented in species tree."
    end #if
    
    @gene_mapping = {}
    @species_numbering = {}
    @spec_node_map = {}
    @gene_node_map = {}
   

  end #initialize

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
    @gene_node_map[@gene_mapping[root]] = root

    if (@gene_mapping[root] == @gene_mapping[nextNodes[0]]) || (@gene_mapping[root] == @gene_mapping[nextNodes[1]])
      root.events = Bio::PhyloXML::Events.new
      root.events.duplications = 1
    else
      root.events = Bio::PhyloXML::Events.new
      root.events.speciations = 1
    end #if

  end #computeMapping

  def _computeMapping(node) 
     
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
      @gene_node_map[gene_mapping[node]] = node
      if (@gene_mapping[node] == @gene_mapping[nextNodes[0]]) || (@gene_mapping[node] == @gene_mapping[nextNodes[1]])
        node.events = Bio::PhyloXML::Events.new
        node.events.duplications = 1
      else
        node.events = Bio::PhyloXML::Events.new
        node.events.speciations = 1
      end #if

    end #if

  end #_computeMapping

  def isRooted?(tree)
    if tree.root == nil:
      return false
    else
      return true
    end #if
    
  end #isRooted?
 
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

  def initializeMap(tree, name_node_map)
    root = tree.root
    index = 1
    map = {}
    map[root] = index
    name_node_map[index] = root
    index += 1
    
    nextNodes = tree.children(root)
    nextNodes.each { |n|
      index, map, name_node_map = _initializeMap(tree, n, index, map, name_node_map)
    }
   
    return map, name_node_map
  end #initializeMap
 
  def _initializeMap(tree, node, index, map, name_node_map)
   
    map[node] = index
    name_node_map[index] = node
    index = index + 1
   
    nextNodes = tree.children(node)
    nextNodes.each{ |n|
      index, map, name_node_map = _initializeMap(tree, n, index, map, name_node_map)
    }
    return index, map, name_node_map

  end #_initializeMap

  def initializeMapping
    
    @species_numbering, @spec_node_map = initializeMap(@species_tree, @spec_node_map)
    index = 0
    @gene_tree.leaves.each{ |n|
      @species_tree.leaves.each{ |s|
        if nodeEqual?(n, s)
          @gene_mapping[n] = @species_numbering[s]
          @gene_node_map[@species_numbering[s]] = n
        end # if
       }
    }
 
  end #initializeMapping

  def nodeEqual?(node1, node2)
  
    if (node1.taxonomies[0].taxonomy_id != nil) && (node2.taxonomies[0].taxonomy_id != nil)
      if (node1.taxonomies[0].taxonomy_id.provider != nil) && (node2.taxonomies[0].taxonomy_id.provider != nil)
        if node1.taxonomies[0].taxonomy_id.value == node2.taxonomies[0].taxonomy_id.value
          return true
        else 
          return false
        end # if
      end #if
  
    elsif (node1.taxonomies[0].code != nil) && (node2.taxonomies[0].code != nil) 
      if node1.taxonomies[0].code == node2.taxonomies[0].code
        return true
      else 
        return false
      end #if
  
    elsif (node1.taxonomies[0].scientific_name != nil) && (node2.taxonomies[0].scientific_name != nil)
      if node1.taxonomies[0].scientific_name == node2.taxonomies[0].scientific_name
        return true
      else 
        return false
      end #if
  
    elsif (node1.taxonomies[0].common_names[0] != nil) && (node2.taxonomies[0].common_names[0] != nil)
      if node1.taxonomies[0].common_names[0] == node2.taxonomies[0].common_names[0]
        return true
      else 
        return false
      end # if
  
    else 
      raise "Nodes must share an identifier to be compared."
    end #if

  end #nodeEqual?

  def writeUpdatedGeneTreePhyloXML(filename)
   
    # Create new phyloxml writer
    writer = Bio::PhyloXML::Writer.new(filename)
   
    # Write tree to the file tree.xml
    writer.write(@gene_tree)

  end #writeUpdatedGeneTreePhyloXML

  def computeSpeciationDuplications

    initializeMapping
    computeMapping

  end #computeSpeciationDuplications

  private :_initializeMap
  private :_computeMapping

end

sdi = SDI.new('species.xml', 'gene.xml')
sdi.computeSpeciationDuplications 

sdi.writeUpdatedGeneTreePhyloXML('test_updated_function.xml')
