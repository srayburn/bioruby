#
# = test/unit/bio/util/phylogeny/SDI/test_sdir.rb - Unit tests for Bio::Algorithm::SDIR
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
require 'bio/util/phylogeny/SDI/SDIR'

module Bio
class TestSDIRData
  SDIR_TEST_DATA = Pathname.new(File.join(BioRubyTestDataPath, 'SDI')).cleanpath.to_s
  def self.gene_xml
    File.join SDIR_TEST_DATA, 'gene.xml'
  end
  def self.species_xml
    File.join SDIR_TEST_DATA, 'species.xml'
  end
  def self.results_xml
    File.join SDIR_TEST_DATA, 'test_results.xml'
  end
  def self.gene_made_up_xml 
    File.join SDIR_TEST_DATA, 'test_G.xml'
  end
  def self.species_made_up_xml
    File.join SDIR_TEST_DATA, 'test_S.xml'
  end
end # class TestSDIData

class TestSDIRClassMethods < Test::Unit::TestCase
   def setup
    @g = Bio::PhyloXML::Parser.open(TestSDIData.gene_xml)
    @g = @g.next_tree
    @s = Bio::PhyloXML::Parser.open(TestSDIData.species_xml)
    @s = @s.next_tree
  end #setup

  def test_new
    sdi = Bio::Algorithm::SDIR.new(@g, @s)
    assert_instance_of(Bio::PhyloXML::Tree, sdi.gene_tree)
    assert_instance_of(Bio::PhyloXML::Tree, sdi.species_tree)
  end #test_new
  
  def test_set_rooted
    sdi = Bio::Algorithm::SDIR.new(@g, @s)
    sdi.set_rooted(sdi.gene_tree, true)
    assert(sdi.gene_tree.root)
    sdi.set_rooted(sdi.species_tree, true)
    assert(sdi.species_tree.root)
  end #test_set_rooted

end
end