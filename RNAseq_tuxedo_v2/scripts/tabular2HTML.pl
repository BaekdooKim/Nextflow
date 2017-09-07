
use strict;

my ($table, $input, $output, @spans) = @ARGV;

# positions that should have dark right borders
#
my %dark;

foreach my $span ( @spans )
{
        if ( $span > 1 )
        {
                # only initialize if there are any spans more than 1

                pop @spans;
                my $i;

                foreach my $span ( @spans )
                {
                $i += $span;
                    $dark{$i - 1} = 1;
                }

                last;
    }
}

open TABLE, ">>$table" or die $!;
open INPUT, "<$input" or die $!;

my $parity;

while ( my $line = <INPUT> )
{
        chomp $line;

        my @vals = split /\t/, $line;

        $parity = ( $parity eq 'odd' ? 'even' : 'odd' );

        print TABLE "<tr class='$parity'>";

    for( my $i = 0; $i < @vals; $i++ )
    {
                my $val = $vals[$i];
                my $class;

                if ( $dark{$i} )
                {
                        $class = ' class="darkRight"';
                }

                $val =~ s/</&lt;/;
                $val =~ s/>/&gt;/;

                if ( $val =~ /https:/ ) {
                    if ($val =~ /genome.ucsc/ ) {
		    	print TABLE "<td$class><a href='$val'>UCSC</a></td>";
                    }		
		    elsif ( $val =~ /ncbi.nlm.nih.gov/ ) {
                        print TABLE "<td$class><a href='$val'>dbSNP</a></td>";
                    }
                    elsif ( $val =~ /genome.ewha/ ) {
                        print TABLE "<td$class><a href='$val'>ECgene</a></td>";
                    }
		    elsif ( $val =~ /genecards.weizmann/ ) {
                        print TABLE "<td$class><a href='$val'>GeneLoc</a></td>";
		    }
		    elsif ( $val =~ /informatics.jax.org/ ) {
                        print TABLE "<td$class><a href='$val'>MGI</a></td>";
                    }
                    elsif ( $val =~ /rgd.mcw.edu/ ) {
                        print TABLE "<td$class><a href='$val'>RatGB</a></td>";
                    }
                    elsif ( $val =~ /genome.jp/ ) {
                        print TABLE "<td$class><a href='$val'>Kegg</a></td>";
                    }
                    elsif ( $val =~ /genecards.org/ ) {
                        print TABLE "<td$class><a href='$val'>GeneCards</a></td>";
                    } 
                    elsif ($val =~ /ensembl.org/ ) {
                        print TABLE "<td$class><a href='$val'>Ensembl</a></td>";
                    } 
                    elsif ($val =~ /string-db.org/ ) {
                        print TABLE "<td$class><a href='$val'>StringDB</a></td>";
                    }
		    else {
                        print TABLE "<td$class><a href='$val'>$val</a></td>";
                    }
		}
		elsif ( $val =~ /http:/ ) {
		    if ( $val =~ /ghr.nlm/ ) {
			print TABLE "<td$class><a href='$val'>GHR</a></td>";
		    }
		    elsif ( $val =~ /genecards.org/ ) {
                        print TABLE "<td$class><a href='$val'>GeneCards</a></td>";
                    }
                    elsif ( $val =~ /ncbi.nlm.nih.gov/ ) {
                        print TABLE "<td$class><a href='$val'>dbSNP</a></td>";
                    }
		    elsif ( $val =~ /drugbank.ca/ ) {
                        print TABLE "<td$class><a href='$val'>DrugBank</a></td>";
                    } 
                    elsif ( $val =~ /genome.jp/ ) {
                        print TABLE "<td$class><a href='$val'>Kegg</a></td>";
                    }
		    elsif ( $val =~ /informatics.jax.org/ ) {
                        print TABLE "<td$class><a href='$val'>MGI</a></td>";
                    }
                    elsif ( $val =~ /pharmgkb.org/ ) {
                        print TABLE "<td$class><a href='$val'>PharmGKB</a></td>";
                    }
                    elsif ($val =~ /genome.ucsc/ ) {
                        print TABLE "<td$class><a href='$val'>UCSC</a></td>";
                    }
                    elsif ($val =~ /genecards.weizmann/ ) {
                        print TABLE "<td$class><a href='$val'>GeneLoc</a></td>";
                    } 
                    elsif ($val =~ /ensembl.org/ ) {
                        print TABLE "<td$class><a href='$val'>Ensembl</a></td>";
                    }
		    elsif ($val =~ /string-db.org/ ) {
                        print TABLE "<td$class><a href='$val'>StringDB</a></td>";
                    }
		    else {
                        print TABLE "<td$class><a href='$val'>$val</a></td>";
                    }
		}
               else {
                    print TABLE "<td$class>$val</td>";
                }  
        }
	 print TABLE "</tr>\n";
}

close INPUT;

print TABLE "</table></body></html>\n";

close TABLE;

`rsync $table $output`;

