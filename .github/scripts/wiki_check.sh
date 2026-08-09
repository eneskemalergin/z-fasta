#!/usr/bin/env bash
set -euo pipefail

wiki_dir=${1:-wiki}

if [[ ! -d "$wiki_dir" ]]; then
  echo "wiki directory not found: $wiki_dir" >&2
  exit 1
fi

perl - "$wiki_dir" <<'PERL'
use strict;
use warnings;
use File::Basename qw(basename);
use File::Spec;

my ($root) = @ARGV;
opendir(my $dh, $root) or die "cannot open $root: $!\n";
my @entries = sort grep { $_ ne q{.} && $_ ne q{..} } readdir($dh);
closedir($dh);

my %pages;
my $failed = 0;

for my $entry (@entries) {
    my $path = File::Spec->catfile($root, $entry);

    if (!-f $path || -l $path) {
        print STDERR "wiki entries must be regular top-level files: $path\n";
        $failed = 1;
        next;
    }
    if ($entry !~ /\.md\z/) {
        print STDERR "wiki source currently permits Markdown pages only: $path\n";
        $failed = 1;
        next;
    }
    if ($entry =~ m{[\\/:*?"<>|]}) {
        print STDERR "wiki filename contains a non-portable character: $path\n";
        $failed = 1;
    }
    if (!-s $path) {
        print STDERR "wiki page is empty: $path\n";
        $failed = 1;
    }

    (my $page = $entry) =~ s/\.md\z//;
    my $folded = lc $page;
    if (exists $pages{$folded}) {
        print STDERR "case-insensitive wiki page collision: $pages{$folded} and $path\n";
        $failed = 1;
    } else {
        $pages{$folded} = $path;
    }
}

for my $required (qw(Home _Sidebar _Footer)) {
    if (!exists $pages{lc $required}) {
        print STDERR "required wiki page is missing: $root/$required.md\n";
        $failed = 1;
    }
}

for my $entry (@entries) {
    next if $entry !~ /\.md\z/;
    my $path = File::Spec->catfile($root, $entry);
    next if !-f $path;

    open(my $fh, q{<}, $path) or die "cannot read $path: $!\n";
    local $/;
    my $text = <$fh>;
    close($fh);

    while ($text =~ /(?<!!)\[[^\]]+\]\(([^)]+)\)/g) {
        my $target = $1;
        $target =~ s/\s+"[^"]*"\z//;
        next if $target =~ m{^(?:https?://|mailto:|#)};

        if ($target =~ /\.md(?:#|\z)/) {
            print STDERR "$path uses a repository-style .md wiki link: $target\n";
            $failed = 1;
        }

        $target =~ s/#.*\z//;
        $target =~ s{^\./}{};
        $target =~ s/\.md\z//;
        next if $target eq q{};

        if (!exists $pages{lc $target}) {
            print STDERR "$path links to a missing wiki page: $target\n";
            $failed = 1;
        }
    }
}

exit 1 if $failed;
print "validated " . scalar(keys %pages) . " wiki pages in $root\n";
PERL
