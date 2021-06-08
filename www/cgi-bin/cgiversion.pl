#!/usr/bin/perl
use CGI;
$page = new CGI;
print $page->header;
print $page->start_html("CGI.pm version information");
print "The CGI.pm version is ",$CGI::VERSION,"\n";
print $page->end_html;
