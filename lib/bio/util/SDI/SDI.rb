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
    @gene_tree = Bio::PhyloXML::Parser.open('G_test.xml')
    @species_tree = Bio::PhyloXML::Parser.open('S_test.xml')
    @gene_tree = @gene_tree.next_tree
    @species_tree = @species_tree.next_tree
    #test rooted
    #need to make exception, with boolean returning functions
    isRooted(@gene_tree)
    isRooted(@species_tree)
    #test subset
    #need to make exception, with boolean returning functions
    isSubset(@species_tree, @gene_tree)
    #map
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

  def preOrder(root, map)
    index = 1
    map[index] = root.taxonomies[0].taxonomy_id.value
      index += 1
      nextNodes = children(root)
      nextNodes.each { |n|
        index = _preOrder(n, index, map)
      }
  end
  def _preOrder(node, index, map)
    map[index] = node.taxonomies[0].taxonomy_id.value
    index += 1
    nextNodes = children(node)
    nextNodes.each{ |n|
      index = _preOrder(n, index, map)
    }
    return index
  end
  private :_preOrder
end

sdi = SDI.new

