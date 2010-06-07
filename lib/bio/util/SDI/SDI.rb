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
    @map = {}
    @map = initializeMapping(@species_tree, @gene_tree)
    puts @map
  end #initialize

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

  def preOrder(tree)
    root = tree.root
    index = 1
    map = {}
    map['root'] = index
    index += 1
    nextNodes = tree.children(root)
    nextNodes.each { |n|
      index = _preOrder(tree, n, index, map)
    }
   
    return map
  end
 
  def _preOrder(tree, node, index, map)
    puts node.taxonomies
    if node.taxonomies[0]:
      map[node.taxonomies[0].taxonomy_id.value] = index
    else
      map['internal'] = index
    end #if
    index += 1
   # puts node
    nextNodes = tree.children(node)
  #  puts nextNodes
    nextNodes.each{ |n|
      index = _preOrder(tree, n, index, map)
    }
    return index
  end

  def initializeMapping(spec_tree, gene_tree)
    species_numbering = {}
    species_numbering = preOrder(spec_tree)
    puts species_numbering
    map = {}
    gene_tree.leaves.each{ |n|
      map[n.taxonomies[0].taxonomy_id.value] = species_numbering[n.taxonomies[0].taxonomy_id.value]
    }
    return map
  end
  
  private :_preOrder

end

sdi = SDI.new

