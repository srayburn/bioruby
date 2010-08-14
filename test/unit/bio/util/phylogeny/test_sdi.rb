#
# = test/unit/bio/util/phylogeny/SDI/test_sdi.rb - Unit tests for Bio::Util::SDI
#
# Copyright::   Copyright (C) 2010
#               Sara Rayburn <sararayburn@gmail.com>
# License::     The Ruby License
#

#  load libraries required
#load 'test/bioruby_test_helper.rb'
require 'pathname'
load Pathname.new(File.join(File.dirname(__FILE__), ['..'] * 4,
			'bioruby_test_helper.rb')).cleanpath.to_s
require 'test/unit'
require 'bio/db/phyloxml/phyloxml_parser'
require 'bio/util/phylogeny/SDI/SDI'

module Bio
class TestSDIData
  SDI_TEST_DATA = Pathname.new(File.join(BioRubyTestDataPath, 'SDI')).cleanpath.to_s
  def self.gene_xml
    File.join SDI_TEST_DATA, 'gene.xml'
  end
  def self.species_xml
    File.join SDI_TEST_DATA, 'species.xml'
  end
  def self.results_xml
    File.join SDI_TEST_DATA, 'test_results.xml'
  end
  def self.gene_made_up_xml 
    File.join SDI_TEST_DATA, 'test_G.xml'
  end
  def self.species_made_up_xml
    File.join SDI_TEST_DATA, 'test_S.xml'
  end
end # class TestSDIData

class TestSDIClassMethods < Test::Unit::TestCase
  
  def setup
    @g = Bio::PhyloXML::Parser.open(TestSDIData.gene_xml)
    @g = @g.next_tree
    @s = Bio::PhyloXML::Parser.open(TestSDIData.species_xml)
    @s = @s.next_tree
  end

  def test_new
    sdi = Bio::Algorithm::SDI.new(@g, @s)
    assert_instance_of(Bio::PhyloXML::Tree, sdi.gene_tree)
    assert_instance_of(Bio::PhyloXML::Tree, sdi.species_tree)
  end
  
  def test_is_rooted
    sdi = Bio::Algorithm::SDI.new(@g, @s)
    assert(sdi.is_rooted?(@g))
    assert(sdi.is_rooted?(@s))
  end
  
  def test_initialize_mapping
    sdi = Bio::Algorithm::SDI.new(@g, @s)
    sdi.initialize_mapping!
    leaves = sdi.species_tree.leaves
     leaves.each { |l|
      max = sdi.species_tree.number_of_nodes + 1
      while l != sdi.species_tree.root 
        map = sdi.species_numbering[l]
        assert( map < max)
        assert_equal(l, sdi.spec_node_map[sdi.species_numbering[l]])
        max = map
        l = sdi.species_tree.parent(l)
      end
      assert(sdi.species_numbering[sdi.species_tree.root] == 1)
    }
    assert_not_nil(sdi.gene_mapping)  
  end

  def test_compute_mappping
    sdi = Bio::Algorithm::SDI.new(@g, @s)
    sdi.initialize_mapping!
    @gene_results = sdi.compute_mapping
    @actual_results = Bio::PhyloXML::Parser.open(TestSDIData.results_xml)
    @actual_results = @actual_results.next_tree
    a_leaves = @actual_results.leaves
    g_leaves = @gene_results.leaves
    assert_equal(@actual_results.root.events.speciations, sdi.gene_tree.root.events.speciations)
    assert_equal(@actual_results.root.events.duplications, sdi.gene_tree.root.events.duplications)
    a_children = @actual_results.children(@actual_results.root)
    t_children = sdi.gene_tree.children(sdi.gene_tree.root)
          
  end

  
  end #class TestSDIClassMethods

 end