## Mitomaster
Source code for Mitomaster http://mitomaster.mitomap.org/

#### Stuff you would need to install this web app should you be so foolish to try
   * [Mitomaster perl module](https://github.com/leipzig/Bio_Mitomaster)
   * CLI Haplogrep.jar - contact hansi.weissensteiner@i-med.ac.at for this
   * [Twitter Bootstrap](http://getbootstrap.com/)
   * [DataTables](https://github.com/DataTables/DataTables)
   * [TableTools](https://github.com/DataTables/TableTools) 
   * [GenomeTools](https://github.com/genometools/genometools) - put in "gt/"
   * Create top level temp directories:  `uploads/` `tmphaplo/` `tmpgt/` `tmp_gff3/`
   * Load mitomap and mitomasterb schemas into postgres database called mito
   * Fill in cgi-bin/Mitoprocess/LocalSettings.pm using sample file
