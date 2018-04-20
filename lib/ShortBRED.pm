
package ShortBRED;

use Exporter;

@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(getAbundanceData expandUniRefIds);

#our ($IdentifyScript, $QuantifyScript, $ParseSSNScript);
#
#$IdentifyScript = "shortbred_identify.py";
#$QuantifyScript = "shortbred_quantify.py";
#$ParseSSNScript = "parse_ssn.py";





sub getAbundanceData {
    my $protFile = shift;
    my $clustFile = shift;
    my $cleanId = shift;
    my $isMerged = shift; # Load a merged results file as opposed to an individual quantify run file

    $cleanId = 1 if not defined $cleanId; # By default we remove the cluster number from the front of the protein name
    $isMerged = 0 if not defined $isMerged;

    my $abd = {metagenomes => [], proteins => {}, clusters => {}};

    if (defined $protFile and -f $protFile) {
        open PROT, $protFile or die "Unable to open protein file $protFile: $!";

        my $header = <PROT>;
        chomp($header);
        my ($feat, @mg) = split(m/\t/, $header);
        push(@{$abd->{metagenomes}}, @mg);

        while (<PROT>) {
            chomp;
            #my ($feature, @mgRes) = split(m/\t/);
            my (@parts) = split(m/\t/);
            my ($clusterNum, $feature, @mgRes);
            if ($isMerged) {
                ($clusterNum, $feature, @mgRes) = @parts;
            } else {
                ($feature, @mgRes) = @parts;
                my $tempId;
                ($clusterNum, $tempId) = split(m/\|/, $feature);
                $feature = $tempId if $cleanId and defined $tempId and $tempId =~ m/^[A-Z0-9]{6,10}$/;
            }
     
            #$feature =~ s/^([^\|]+)\|// if $cleanId;
            for (my $i = $#mgRes; $i < $#mg; $i++) { # Ensure that there are the same amount of results as metagenome headers
                push(@mgRes, 0);
            }
     
            #push(@{$abd->{proteins}->{$feature}}, @mgRes);
            for (my $i = 0; $i <= $#mg; $i++) {
                my $mgId = $mg[$i];
                $abd->{proteins}->{$feature}->{$mgId} = $mgRes[$i];
            }
        }

        close PROT;
    }

    if (defined $clustFile and -f $clustFile) {
        open CLUST, $clustFile or die "Unable to open cluster file $clustFile: $!";

        my $header = <CLUST>;
        chomp($header);
        my ($feat, @mg) = split(m/\t/, $header);
        push(@{$abd->{metagenomes}}, @mg) if not scalar @{$abd->{metagenomes}};

        while (<CLUST>) {
            chomp;
            my ($feature, @mgRes) = split(m/\t/);
            for (my $i = $#mgRes; $i < $#mg; $i++) { # Ensure that there are the same amount of results as metagenome headers
                push(@mgRes, 0);
            }
            #push(@{$abd->{clusters}->{$feature}}, @mgRes);
            for (my $i = 0; $i <= $#mg; $i++) {
                my $mgId = $mg[$i];
                $abd->{clusters}->{$feature}->{$mgId} = $mgRes[$i];
            }
        }

        close CLUST;
    }

    return $abd;
}


sub expandUniRefIds {
    my $nodeId = shift;
    my $xmlNode = shift;
    my $efiAnnoUtil = shift;


    my @nodes;

    my @annotations = $xmlNode->findnodes('./*');

    foreach my $annotation (@annotations) {
        my $attrName = $annotation->getAttribute('name');
        if ($efiAnnoUtil->is_expandable_attr($attrName)) {
            my @accessionlists = $annotation->findnodes('./*');
            foreach my $accessionlist (@accessionlists) {
                #make sure all accessions within the node are included in the gnn network
                my $attrAcc = $accessionlist->getAttribute('value');
                #print "Expanded $nodeId into $attrAcc\n";
                push @nodes, $attrAcc if $nodeId ne $attrAcc;
            }
        }
    }

    return @nodes;
}


1;

