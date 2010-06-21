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
    @gene_tree = Bio::PhyloXML::Parser.open('species.xml')
    @species_tree = Bio::PhyloXML::Parser.open('gene.xml')
    @gene_tree = @gene_tree.next_tree
    @species_tree = @species_tree.next_tree
    
    # test rooted
    # need to make exception, with boolean returning functions
    isRooted(@gene_tree)
    isRooted(@species_tree)
    
    # test subset
    # need to make exception, with boolean returning functions
    isSubset?(@species_tree, @gene_tree)
    
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

  def postOrder #(g_tree, s_tree, map, spec_nodes, gene_nodes)
    puts @gene_tree.class
    puts @species_tree.class
    root = @gene_tree.root
    nextNodes = @gene_tree.children(root)
    nextNodes.each { |n|
      index = _postOrder(n) 
    }
    # this is how forester does it…
    a = @gene_mapping[nextNodes[0]]
    b = @gene_mapping[nextNodes[1]]
    
    @gene_mapping[root] = @species_numbering[@species_tree.lowest_common_ancestor(@spec_node_map[a], @spec_node_map[b])]
    #puts @gene_mapping[root]
    @gene_node_map[@gene_mapping[root]] = root
    #puts @gene_node_map[@gene_mapping[root]]
    if (@gene_mapping[root] == @gene_mapping[nextNodes[0]]) || (@gene_mapping[root] == @gene_mapping[nextNodes[1]])
      puts "duplication"        
      root.events = Bio::PhyloXML::Events.new
      root.events.duplications = 1
      puts root.events.speciations
    else
      puts "speciation"
      root.events = Bio::PhyloXML::Events.new
      root.events.speciations = 1
      puts root.events.speciations
    end #if
     
   

  end #postorder

  def _postOrder(node) #tree, s_tree, node,map, spec_nodes, gene_nodes)
     
    nextNodes = @gene_tree.children(node)
    nextNodes.each { |n|
      index = _postOrder(n) # tree, s_tree, n, map, spec_nodes, gene_nodes)
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
        puts "duplication"
        node.events = Bio::PhyloXML::Events.new
        node.events.duplications = 1
        puts node.events.duplications
      else
        puts "speciation"
        node.events = Bio::PhyloXML::Events.new
        node.events.speciations = 1
        puts node.events.speciations
      end #if

    end #if
    #puts index
    #index += 1
    #return index

  end #_postOrder

  def isRooted(tree)
    if tree.root == nil:
      exit
    end #if
    
  end #isRooted
 
  def isSubset?(superset, subset)
    superset_leaves = superset.leaves

    subset_leaves = subset.leaves
    
    subset_leaves.each { |node|
      puts "NODE OUTER:"
      puts node.inspect
      inner = false
      superset_leaves.each{ |node1|
        puts "NODE INNER"
        puts node1.inspect
        inner = false
        if nodeEqual?(node, node1)
          inner = true
          break
        end #if
      }
      if !inner
        puts "Error"
        exit
      end #if
    }
          #puts species.find { |spec| spec == key }
	#if node.taxonomies[0].taxonomy_id.value in species:
	#	puts "ok"
	#end
    
    end #isSubset
 
  def getNodeName(node, index)
    #bio_node = node.to_biotreenode
    if node.taxonomies[0] == nil:
      name = index
    else
      name = node.taxonomies[0].code
    end #if
  end #getNodeName

  def preOrder(tree, name_node_map)
    root = tree.root
    #puts root.class
    index = 1
    #name = getNodeName(root, index)
    #key = root.taxonomies[0].taxonomy_id 
    map = {}
    map[root] = index
    name_node_map[index] = root
   # puts index.class
    #puts 1.class
    index += 1
    
    nextNodes = tree.children(root)
    nextNodes.each { |n|
      index, map, name_node_map = _preOrder(tree, n, index, map, name_node_map)
    }
   
    return map, name_node_map
  end
 
  def _preOrder(tree, node, index, map, name_node_map)
   
    map[node] = index
    name_node_map[index] = node
    index = index + 1
   
    nextNodes = tree.children(node)
    nextNodes.each{ |n|
      index, map, name_node_map = _preOrder(tree, n, index, map, name_node_map)
    }
    return index, map, name_node_map
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
      species_tree.leaves.each{ |s|
        if nodeEqual?(n, s)
          gene_mapping[n] = species_numbering[s]
          gene_node_map[species_numbering[s]] = n
          puts "equal"
        end # if
       }
       if (gene_mapping[n] == nil)
         puts "Error mapping"
         exit
       end # if
    }
 
    return gene_mapping, species_numbering, spec_node_map, gene_node_map
  end

  def nodeEqual?(node1, node2)
    puts node1.taxonomies[0]
    puts node2.taxonomies[0]
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
      puts "Error"
      exit
    end #if

  end #nodeEqual?

  def writeUpdatedGeneTreePhyloXML(filename)
    # Create new phyloxml writer
    writer = Bio::PhyloXML::Writer.new(filename)
    # Write tree to the file tree.xml
    writer.write(@gene_tree)
  end #writeUpdatedGeneTreePhyloXML

  private :_preOrder
  private :_postOrder

end

sdi = SDI.new
sdi.postOrder #(sdi.gene_tree, sdi.species_tree, sdi.gene_mapping, sdi.spec_node_map, sdi.gene_node_map)
root = sdi.gene_tree.root
puts root.class
puts root.events
puts root.events.speciations
sdi.writeUpdatedGeneTreePhyloXML('test_updated_function.xml')
