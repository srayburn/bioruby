#!/usr/bin/env ruby
# = util/SDI/sdi_binary.rb - Binary SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#

#  load libraries required
$: << '/Users/sarara/bioruby/lib'
require 'bio/db/phyloxml/phyloxml_parser'

module Bio
  module SDI

  #  def load
      @G = Bio::PhyloXML::Parser.open('G_test.xml')
      @S = Bio::PhyloXML::Parser.open('S_test.xml')
   # end

    #def isRooted
    @S_tree = @S.next_tree
    if @S_tree.root == nil:
	exit
    end
    @G_tree = @G.next_tree
    if @G_tree.root == nil:
        exit
    end
   # end

   # def isSubset
      @S_tree = @S.next_tree
      S_leaves = @S_tree.leaves
      species = []
      S_leaves.each { |node|
       species[species.length] = node.taxonomies[0].taxonomy_id.value
      }


      @G_tree = @G.next_tree
      G_leaves = @G_tree.leaves
      G_leaves.each { |node|
	key = node.taxonomies[0].taxonomy_id.value
	puts species.find { |spec| spec == key }
	#if node.taxonomies[0].taxonomy_id.value in species:
	#	puts "ok"
	#end
      }

    def preOrder(root, map)
      index = 1
      map[index++] = root.taxonomies[0].taxonomy_id.value
      nextNodes = children(root)
      nextNodes.each { |n|
        index = _preOrder(n, index, map)
    end
    def _preOrder(node, index, map)
      map[index++] = node.taxonomies[0].taxonomy_id.value
      nextNodes = children(node)
      nextNodes.each{ |n|
        index = _preOrder(n, index, map)
      return index
    end
    private :_preOrder(node)
    end
end
end