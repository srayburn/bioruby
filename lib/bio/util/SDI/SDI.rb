# = util/SDI/sdi_binary.rb - Binary SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#

#  load libraries required
#$: << '/Users/sarara/bioruby/lib'
require 'bio/db/phyloxml/phyloxml_parser'


class SDI

  attr_accessor :gene_tree
  attr_accessor :species_tree
  attr_accessor :gene_mapping 
  attr_accessor :species_numbering
  attr_accessor :spec_node_map
  attr_accessor :gene_node_map

  def initialize
    # load phyloxml trees (assumption: each file contains only 1 tree)
    @gene_tree = Bio::PhyloXML::Parser.open('G_test.xml')
    @species_tree = Bio::PhyloXML::Parser.open('S_test.xml')
    @gene_tree = @gene_tree.next_tree
    @species_tree = @species_tree.next_tree
    
    # test rooted
    # need to make exception, with boolean returning functions
    isRooted(@gene_tree)
    isRooted(@species_tree)
    
    # test subset
    # need to make exception, with boolean returning functions
    isSubset(@species_tree, @gene_tree)
    
    # map
    @gene_mapping = {}
    @species_numbering = {}
    @spec_node_map = {}
    @gene_node_map = {}
    @gene_mapping, @species_numbering, @spec_node_map, @gene_node_map = initializeMapping(@species_tree, @gene_tree)
    puts @gene_mapping
    puts @spec_node_map
    puts @gene_node_map
  end #initialize

  def postOrder(tree, map, spec_nodes)
    puts tree.class
    root = tree.root
    nextNodes = tree.children(root)
    nextNodes.each { |n|
      index = _postOrder(tree, n, index, map, spec_nodes)
    }
    # this is how forester does it…
    a = map[getNodeName(nextNodes[0], index)]
    b = map[getNodeName(nextNodes[1])]
    
    map[root] = LCA(spec_nodes[a], spec_nodes[b])
    puts map[root]
    
     
   

  end #postorder

  def _postOrder(tree, node, index, map, spec_nodes, gene_nodes)
     
    nextNodes = tree.children(node)
    nextNodes.each { |n|
      index = _postOrder(tree, n, map, spec_nodes)
    }
    if nextNodes[0] != nil:
      a = map[getNodeName(nextNodes[0], index)]
      if nextNodes[1] != nil:
        b = map[getNodeName(nextNodes[1], index)]
      else 
        b = map[getNodeName(nextNodes[0], index)]
      end #if  
    
      map[node.getNodeName] = LCA(spec_nodes[a], spec_nodes[b])
      puts map[getNodeName(node), index]
    end #if
    index += 1
    return index

  end #_postOrder

  def isRooted(tree)
    if tree.root == nil:
      exit
    end #if
    
  end #isRooted
 
  def isSubset(superset, subset)
    superset_leaves = superset.leaves
    species = []
    superset_leaves.each { |node|
      species[species.length] = node.taxonomies[0].taxonomy_id.value
    }

    subset_leaves = subset.leaves
    subset_leaves.each { |node|
      key = node.taxonomies[0].taxonomy_id.value
      puts species.find { |spec| spec == key }
	#if node.taxonomies[0].taxonomy_id.value in species:
	#	puts "ok"
	#end
    }
    end #isSubset
 
  def getNodeName(node, index)
    #bio_node = node.to_biotreenode
    if node.taxonomies[0] == nil:
      name = index
    else
      name = node.taxonomies[0].taxonomy_id.value
    end #if
  end #getNodeName

  def preOrder(tree, name_node_map)
    root = tree.root
    #puts root.class
    index = 1
    name = getNodeName(root, index)

    map = {}
    map[name] = index
    name_node_map[index] = root
   # puts index.class
    #puts 1.class
    index += 1
    
    nextNodes = tree.children(root)
    nextNodes.each { |n|
      index, name_node_map = _preOrder(tree, n, index, map, name_node_map)
    }
   
    return map, name_node_map
  end
 
  def _preOrder(tree, node, index, map, name_node_map)
    #puts node.taxonomies
    name = getNodeName(node, index)
    map[name] = index
    name_node_map[index] = node
    #puts index
    #puts 1.class
     index = index + 1
   
   # puts node
    nextNodes = tree.children(node)
  #  puts nextNodes
    nextNodes.each{ |n|
      index, name_node_map = _preOrder(tree, n, index, map, name_node_map)
    }
    return index, name_node_map
  end

  def initializeMapping(spec_tree, gene_tree)
    species_numbering = {}
    spec_node_map = {}
    species_numbering, spec_node_map = preOrder(spec_tree, spec_node_map)
    puts species_numbering
    gene_mapping = {}
    gene_node_map = {}
    index = 0
    gene_tree.leaves.each{ |n|
      name = getNodeName(n, index)
      index += 1
      gene_mapping[name] = species_numbering[name]
      gene_node_map[gene_mapping[name]] = n
    }
 
    return gene_mapping, species_numbering, spec_node_map, gene_node_map
  end

  private :_preOrder
  private :_postOrder

end

sdi = SDI.new
puts 
sdi.postOrder(sdi.gene_tree, sdi.gene_mapping, sdi.spec_node_map)
