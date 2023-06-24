
use File::Spec;
sub convert_amino_acid_sequence{
    
    my $amino_acid_sequence;
    my $markers;
    my $offset;
    my $id;
    my $separator;
    my $binaries;
    my $positions;
    my $strict;
    my $verbose;
    my $self;
    my $ref = ref $_[0];
    my %parameters;

    if( $ref->isa( __PACKAGE__ ) ){
        $self = shift;
        $strict = $self->{"Strict"};
        $offset = $self->{"Offset"};
        $verbose = $self->{"Verbose"};
        $markers = $self->{"Markers"};
    }

    %parameters = @_;
    
    $amino_acid_sequence = $parameters{"AminoAcidSequence"} if $parameters{"AminoAcidSequence"};
    $markers = $parameters{"Markers"} if $parameters{"Markers"};
    $offset = 0;
    $id = $parameters{"ID"} if $parameters{"ID"};
    $binaries = $parameters{"Binaries"};
    $positions = $parameters{"Positions"};
    $strict = $parameters{"Strict"} if $parameters{"Strict"};
    $verbose = $parameters{"Verbose"} if $parameters{"Verbose"};
        
    if( !$markers ){
        $self->parse_markers_file();
        $markers = $self->{"Markers"};
    }
    
    for my $i ( 0 .. $#{$markers} ){
        
        my $amino_acid;

        if( ref $amino_acid_sequence eq 'HASH'){
            $amino_acid = $amino_acid_sequence->{  ( $markers->[$i]->{"Position"} + $offset ) };
            $amino_acid = $codon_to_amino_acid->{$amino_acid};
        }
        else{
            $amino_acid = substr( $amino_acid_sequence, $markers->[$i]->{"Position"} + $offset, 1 );
        }
        if( !$amino_acid ){
            die "There is no amino acid at " . $markers->[$i]->{"Position"} + $offset . " in sequence $id.\n" if $strict;
            push @$binaries, 0 if $binaries;
            next;
        }
        my $Needle = quotemeta( $amino_acid );
        if( grep /$Needle/, @{$markers->[$i]->{"AminoAcids"}} ){
            push @$binaries, 1 if $binaries;
            push @$positions, $i + 1 if $positions;
        }
        else{
            push @$binaries, 0 if $binaries;
        }
    }
    if( scalar @$positions == 0 and $positions ){
        push @$positions, 0;
    }
    return 1;
}
sub load_sequences{
    
    my $self = shift;
    my %parameters = @_;
    setup_variable( \my $binaries_directory, $parameters{Binaries_Directory}, $self->{Binaries_Directory}, "binaries/" ) and do{
  	say STDOUT "Error: Could not figure out which directory contains sminpd binaries";
        return -1;
    };
    
    return -1 if setup_variable( \my $verbose, $parameters{Verbose}, $self->{Verbose}, 0 );
    return -1 if setup_log( \my $log, $parameters{Log}, \*STDERR );
    return -1 if setup_variable( \my $gap_character, $parameters{Gap_Character}, $self->{Gap_Character}, '-' );
    return -1 if setup_variable( \my $dna_directory, $parameters{DNA_Directory}, $self->{DNA_Directory} );
    return -1 if setup_variable( \my $format, $parameters{DNA_Input_Format}, $self->{DNA_Input_Format}, 'FASTA' );
    
    my $binaries = {};
    my $sequences = {};
    my $positions = {};
    
    $self->parse_markers_file( Verbose => $verbose, Log => $log ) if !$self->{Markers};
    my $markers = $self->{Markers};
    return -1 if !$markers;
    
    if( -d File::Spec->rel2abs( $dna_directory ) ){
        $dna_directory = File::Spec->rel2abs( $dna_directory );
        opendir DNA_Directory, $dna_directory or do{ print $log "Error: Couldn't open \"$dna_directory\" : $!\n"; return -1; };
        $dna_directory = [ map { File::Spec->catfile( $dna_directory, $_ ) } readdir DNA_Directory ];
        closedir DNA_Directory;
    }
    else{
        #~ $dna_directory = [ Tkx::SplitList($dna_directory) ]; ## Can't figure out WHAT this was doing...
    }
    $dna_directory = [ map { File::Spec->rel2abs($_); } @$dna_directory ];
    $dna_directory = [ grep { !( -d $_ )} @$dna_directory ];
    
    say $log "Loading DNA sequences ..." if $verbose;
    
    for my $dna_file ( @$dna_directory ) {
        next if -d $dna_file;
        my $basename = basename($dna_file);
        $basename =~ /^(\d+)\./;
        my $patient = $1;
        say $log "\t$dna_file" if $verbose;
        my %sequences = ();
        my $seq_in = Bio::SeqIO->new(
            -format => $format,
            -file   => $dna_file,
        );
        
        while ( my $seq = $seq_in->next_seq() ) {
            if ( $seq->alphabet ne 'dna' ) {
                print $log "Error: Found non-DNA sequence in file \"$dna_file\"\n";
                next;
            }
            my $display_id = $seq->display_id();
            my $sequence = $seq->seq;
            for my $i ( 0 .. $#{$markers} ){
                my $nucleotide_position = $markers->[$i]->{Position};
                my $amino_acid_position = $nucleotide_position / 3;
                my @codon;
                for my $j ( 0 .. 2 ){
                    my $nucleotide = substr $sequence, $nucleotide_position + $j, 1;
                    push @codon, $nucleotide ? $nucleotide : $gap_character;
                }
                $sequences->{$patient}{$display_id}{DNA}{$nucleotide_position} = join '', @codon;
                $sequences->{$patient}{$display_id}{Protein}{ $amino_acid_position } = $codon_to_amino_acid->{join '', @codon};
            }
            my $id = "$patient.$display_id";
            $binaries->{$patient}{$display_id} = [];
            $positions->{$patient}{$display_id} = [];
            $self->convert_amino_acid_sequence(
                AminoAcidSequence => $sequences->{$patient}{$display_id}{DNA},
                ID => $id,
                Binaries => $binaries->{$patient}{$display_id},
                Positions => $positions->{$patient}{$display_id},
            );
	    #say join ',', ($id, join('', @{$binaries->{$patient}{$display_id}}), join('', @{$positions->{$patient}{$display_id}})) if $verbose > 2;
        }
    }

    $self->{Binaries} = $binaries;
    $self->{Positions} = $positions;
    $self->{Sequences} = $sequences;

    say $log "Done loading sequences" if $verbose;
    
    undef $log;
    return 0;
    
}

sub setup_variable{
    
    my $self = shift;
    unless( ref( $self ) eq __PACKAGE__ || $self eq __PACKAGE__ ){
        unshift @_, $self;
    }
    
    my $ref = shift or return -1;
    
    for my $arg ( @_ ){
        if( defined $arg ){
            $$ref = $arg;
            return 0;
        }
    }
    
    return -1;
}   
sub load_transitions{
    
	#my $self = shift;
    my %parameters = @_;
    while ( ($k,$v) = each %parameters ) {
    print "$k => $v\n";
}
print($parameters{Transition_Directory});
    my $transitions = {};
    
    setup_log( \my $log, $parameters{Log}, $self->{Log}, \*STDERR ) and do{
        say STDOUT "Error: Could not set up a log file for writing";
        return -1;
    };
    #setup_variable( \my $transitions_directory, $parameters{Transition_Directory}, $self->{Transitions_Directory}, File::Spec->catfile( $self->{SlidingMinPD_Directory}, "transitions" ) ) and do{
    setup_variable( \my $transitions_directory, $parameters{Transition_Directory}, $self->{Transition_Directory}, "transitions3/" ) and do{
  	say STDOUT "Error: Could not figure out which directory contains sminpd transitions";
        return -1;
    };
    say STDOUT "Loading transitions ...";
    
    opendir TRANSITIONS_DIRECTORY, $transitions_directory or do{
        say STDOUT "Error: Couldn't open '$transitions_directory': $!";
        return -1;
    };
    say STDOUT "SSS.";
    while( my $base = readdir TRANSITIONS_DIRECTORY ){
	say STDOUT $base;
        next if $base !~ /^(\d+)\./;
	say STDOUT "BYE";
        my $patient_id = $1;
	say STDOUT $patient_id;
	say STDOUT "HERE";
        my $file = File::Spec->catdir( $transitions_directory, $base );
	say STDOUT "HERE";
        next if -d $file;
	say STDOUT "HERE";
        open FILE, $file or return $!;
	say STDOUT "HERE";
        while( my $line = <FILE> ){
	    say STDOUT "MMM";
            chomp $line;
            next if $line !~ /^\d+/;
            my( $descendant, $ancestor, $distance, $bootstrap ) = split / /, $line;
            $transitions->{$patient_id}{$ancestor}{$descendant}{Bootstrap} = $bootstrap;
            say STDOUT join ',', ($patient_id, $ancestor, $descendant, $bootstrap );
        }
        close FILE;
    }
    close TRANSITIONS_DIRECTORY;

    #$self->{Transitions} = $transitions;
    
    say STDOUT "Done loading transitions.";
    
    return $transitions;#0;
    
}
sub setup_log{
    # Returns a FILEHANDLE
    my $self = shift;
    unless( ref( $self ) eq __PACKAGE__ || $self eq __PACKAGE__ ){
        unshift @_, $self;
    }
    my $old_output = shift or return -1;
    my $new_output = "";
    for my $arg ( @_ ){
        if( $arg ){
            $new_output = $arg;
            last;
        }
    }
    return -1 if !$new_output;

    if( ref( $new_output ) eq 'GLOB' ){
        $$old_output = $new_output;
    }
    else{
        open $$old_output, ">>" . $new_output or die $!;
    }
    
    return 0;
    
}
sub setup_output{
    my $self = shift;
    unless( ref( $self ) eq __PACKAGE__ || $self eq __PACKAGE__ ){
        unshift @_, $self;
    }
    my $old_output = shift or return -1;
    my $new_output = "";
    for my $arg ( @_ ){
        if( $arg ){
            $new_output = $arg;
            last;
        }
    }
    return -1 if !$new_output;

    if( ref( $new_output ) eq 'GLOB' ){
        $$old_output = $new_output;
    }
    else{
        open $$old_output, ">" . $new_output or die $!;
    }
    
    return 0;
    
}

sub calculate_summary{
    
	#my $self = shift;
    my %parameters = @_;

    setup_log( \my $log, $parameters{Log}, $self->{Log}, \*STDERR );
    setup_variable( \my $verbose, $parameters{Verbose}, $self->{Verbose}, 0 );
    setup_variable( \my $bootstrap_threshold, $parameters{Bootstrap_Threshold}, $self->{Bootstrap_Threshold}, 0 );
    setup_variable( \my $aggregate, $parameters{Aggregate}, $self->{Aggregate}, 1 );
    say STDERR "MADE IT";   
    my $transitions;
    my $binaries;
    my $summary = {};
    my $transitions = load_transitions( Verbose => $verbose, Log => $log, Transition_Directory => $parameters{Transition_Directory} ); 
    #if !$self->{Transitions};
    #my $transitions = $self->{Transitions};
    if( !$transitions ){
        say STDOUT "Error: Couldn't find transitions."; 
        return -1;
    }
    #load_sequences( Verbose => $verbose, Log => $log ) if !$self->{Binaries};
    my $binaries = load_sequences(Verbose => $verbose, Log => $log, Binaries_Directory => $parameters{Binaries_Directory} ) if !$self->{Binaries}; 
    
    #$self->{Binaries};
    if( !$binaries ){
        say STDOUT "Error: Couldn't find binaries." ;
        return -1;
    }
    
    
    say STDOUT "Calculating summary ...\n" if $verbose;
    
    for my $patient ( sort keys %$transitions ){
        my %aggregates;
        for my $ancestor ( sort keys %{$transitions->{$patient}} ){
            for my $descendant ( sort keys %{$transitions->{$patient}{$ancestor}} ){
                my $bootstrap = $transitions->{$patient}{$ancestor}{$descendant}{Bootstrap};
                if( $bootstrap < $bootstrap_threshold ){
                    say $log "skipping $patient.$ancestor because $bootstrap < $bootstrap_threshold ($verbose)" if $verbose > 1;
                    next;
                }
                else{
                    say $log "not skipping $patient.$ancestor because $bootstrap >= $bootstrap_threshold ($verbose)" if $verbose > 1;
                }
                my $ancestor_binary = join '', @{$binaries->{$patient}{$ancestor}};
                my $descendant_binary = join '', @{$binaries->{$patient}{$descendant}};
                my $ancestor_time = (split /\./, $ancestor)[0];
                my $descendant_time = (split /\./, $descendant)[0];
                my $time_difference = $descendant_time - $ancestor_time;
                my @row_values = (
                    $descendant_binary,
                    $ancestor_binary,
                    $time_difference,
                    $descendant_time,
                    $ancestor_time,
                    $patient,
                    $ancestor,
                );
                push @row_values, $descendant if !$aggregate;
                $aggregates{ join ',', @row_values }++;
            }
        }
        for my $row ( sort keys %aggregates ){
            my( $descendant_binary, $ancestor_binary, $time ) = split /,/, $row;
            $summary->{$ancestor_binary}{$descendant_binary}{Count}++;
            $summary->{$ancestor_binary}{$descendant_binary}{Time} += $time;
        }
    }
    
    $self->{Summary} = $summary;
    
    print $log "Done calculating summary ...\n" if $verbose;
    
    return $summary;
    
}

sub write_summary{
    
	#my $self = shift;
    my %parameters = @_;
    while ( ($k,$v) = each %parameters ) {
    print "$k => $v\n";
}

    setup_log( \my $log, $parameters{Log}, $self->{Log}, \*STDERR );
    setup_output( \my $output, $parameters{Output}, "summary.csv", \*STDOUT );
    setup_variable( \my $headers, $parameters{Headers}, $self->{Headers}, 1 );
    setup_variable( \my $verbose, $parameters{Verbose}, $self->{Verbose}, 0 );
    setup_variable( \my $bootstrap_threshold, $parameters{Bootstrap_Threshold}, $self->{Bootstrap_Threshold}, 0 );
    setup_variable( \my $aggregate, $parameters{Aggregate}, $self->{Aggregate}, 1 );
    
    setup_variable( \my $summary, $parameters{Summary}, $self->{Summary} );

    if( !$summary ){
	    #calculate_summary() if !$self->{Summary};
        $summary = calculate_summary(Verbose => $verbose, Log => $log, Transition_Directory => $parameters{Transition_Directory});
	#$self->{Summary} if $self->{Summary};
    }
    if (!$summary) {
       say STDOUT "NO SUMM"; return -1;
    }
    #return -1 if !$summary;
    
    print $log "Writing summary ...\n" if $verbose;

    my @header_values = (
        'descendant binary',
        'ancestor binary',
        'time',
        'count',
    );
    print $output join( ',', @header_values ) . "\n" if $headers;
    
    for my $ancestor_binary ( sort keys %$summary ){
        for my $descendant_binary ( sort keys %{$summary->{$ancestor_binary}} ){
            my $count = $summary->{$ancestor_binary}{$descendant_binary}{Count};
            my $time = $summary->{$ancestor_binary}{$descendant_binary}{Time};
            my @row_values = (
                $descendant_binary,
                $ancestor_binary,
                $time,
                $count,
            );
            print $output join( ',', @row_values ) . "\n";
        }
    }
    
    print $log "Done writing summary\n" if $verbose;
    
    undef $log;
    undef $output;
    
    return 0;
    
}

sub input {
   $input_dir = @_[0];
}

sub run {
}

sub output {


write_summary(
            Output => @_[0] ,
            Transition_Directory => File::Spec->catfile($input_dir, "transitions"),
            Binaries_Directory => File::Spec->catfile($input_dir, "binaries")
            #Aggregates_File_Output => "aggregates.csv"
        ) and do {
            $self->{errorMessage} = "The 'transitions' report was not generated successfully. Check the log file.";
        };

}

