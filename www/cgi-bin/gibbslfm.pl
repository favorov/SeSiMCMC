#!/usr/bin/perl
#
# $Id: gibbslfm.pl 1720 2012-05-03 17:30:42Z favorov $
#
# The script is cgi support for SeSiMCMC program, (c) A.Favorov 2001-2
#
# It starts with different parameters a set of times.
# 
# action=version prints version info by starting $sampler_name --version
# 
# action=cat tests whether cats the file refereed by tag source=filename
# is html or not, cats while convertting to HTML if necessary.
# Non-html is formatted to 80 chars per line if format80=yes is set in form. 
#
# action= << Recall    
#	creates a redirection to rindex
#
# action=rindex	
# creates a "recalling index page" for given id (id=XXXXXX)
# which has fields: result, FastA, config, errorlog with the
# corresponding calls to action=cat
#
# action= index    
# creates an "index page" for given id (id=XXXXXX)
# which has fields: result, FastA, config, errorlog with the
# corresponding calls to action=cat
#
# action=form
# action=Avanced\ \\\/
# action=Simple\ \/\\
# action=Completion\ example
# action='Footprint search as in DMMPMM '
# creates the main form with paramters, which names and values 
# mainly equals to the 
# parameter names and values of the configuration file. 
# The formtype affects whether the advanced parameters widget is
# hidden or not.
# Example completes the fields with the arginine example, formtype=simple

# action='Run Sampler>>'
# we generate a temp dir name gibbslfm.XXXXXX (XXXXXX is id given to the task)
# then we put FastA.dnc file with fasta data there, config with configuration 
# and, a html file with "wait for results..." string, the task ID and
# refresh=60sec command in meta=string to result
#
# Then, we send index page (action=index) for that id and 
# close the html page and stdout to break the connection.
#
# Then, we starts 
# "$sampler_name --config $dir/config --output-tables --html --cgi --id $ID --ip $ENV{REMOTE_ADDR} --log $WORKDIR/logfile < $dir/FastA.seq > $dir/result"
#
#
#
#

#Work dir
my $WORKDIR  = "./work";
my $sampler_name="./SeSiMCMC";
my $SeSiMCMC_home="http://favorov.bioinfolab.net/SeSiMCMC";
use strict;
use CGI;
# for waitpid()
use POSIX ":sys_wait_h";

my $page=new CGI;

# stop buffering output
$| = 1;

$page->param(-name=>'action',-value=>'form')
	if (!defined $page->param('action') ||
			$page->param('action') eq '');

for ($page->param('action')) 
{
	/version/  				        && do {&print_version($page);last;};
	/form/			         && do {&form($page,"simple");last;};
	/Advanced\ \\\//      && do {&form($page,"advanced");last;};
	/Simple\ \/\\/       && do {&form($page,"simple");last;};
	/Completion\ example/       && do {&form($page,"example");last;};
	/Footprint search as in DMMPMM /   && do {&form($page,"dmmpmm");last;};
	/Run\ sampler\ >>/     && do {&run_sampler($page);last};
	/cat/            && do {&cat($page);last};
	/<<\ Recall/          && do {&makerindexredir($page);last;};
	/rindex/          && do {&makeindex($page,1);last;};
	/index/          && do {&makeindex($page,0);last;};
	
	print $page->header,
				$page->start_html(-title=>"Error. Unknown action",-bgcolor=>'white'),
				$page->h2("Error. Unknown action"),
				$page->end_html;
}

	print $page->header,
				$page->start_html(-title=>"Error. Empty action",-bgcolor=>'white'),
				$page->h2("Error. Empty action"),
				$page->end_html unless $page->param('action');


sub form
{
	my $default_worst_to_retrieve=1;
	my ($Script)=<<ENDSCRIPT;
	function SelectIfDefault(o_name)
	{
		var text_obj=document.MainInputForm.elements[o_name]
		if (text_obj.value=='default') text_obj.select()
	}

	function LeaveIntField(o_name)
	{
		var text_obj=document.MainInputForm.elements[o_name]
		var text=text_obj.value
		if (isNaN(parseInt(text,10))) text='default'
		else text=parseInt(text,10).toString(10)
		text_obj.value=text
	}

	function LeaveFloatField(o_name)
	{
		var text_obj=document.MainInputForm.elements[o_name]
		var text=text_obj.value
		if (isNaN(parseFloat(text))) text='default'
		else text=parseFloat(text).toString()
		text_obj.value=text
	}


	function LeaveProbabilityField(o_name)
	{
		var text_obj=document.MainInputForm.elements[o_name]
		var text=text_obj.value
		var val
		if (isNaN(parseFloat(text))) text='default'
		else {
			val=parseFloat(text)
			if (val >= 0. && val <= 1.) text=val.toString()
			elseif (val == 0 || val == 1) text=val.toString()
			else text='default'
		}
		text_obj.value=text
	}

	
	function SyncBackgrGrey()
	{
		var chk=document.MainInputForm.elements['manual_background']
		if (chk.checked)
		{
			document.MainInputForm.elements['background_A'].value='100'
			document.MainInputForm.elements['background_T'].value='100'
			document.MainInputForm.elements['background_G'].value='100'
			document.MainInputForm.elements['background_C'].value='100'
		}
		else 
		{
			document.MainInputForm.elements['background_A'].value='auto'
			document.MainInputForm.elements['background_T'].value='auto'
			document.MainInputForm.elements['background_G'].value='auto'
			document.MainInputForm.elements['background_C'].value='auto'
		}
			
	}
	
	function LeaveBackgroundField(o_name)
	{
		var chk=document.MainInputForm.elements['manual_background']
		var text_obj=document.MainInputForm.elements[o_name]
		var text=text_obj.value
		if (chk.checked)
		{
			if (isNaN(parseFloat(text))) text='100'
			else text=parseFloat(text).toString()
		}
		else text='auto'
		text_obj.value=text
	}

ENDSCRIPT

	my ($form,$formtype)=@_;
	
	my ($url)=$form->url();
	$url=~s/\?.*//;

	my $is_example=(defined $formtype && ($formtype eq 'example'));
	
	$formtype='simple' if (!defined $formtype || (($formtype ne 'advanced') && ($formtype ne 'dmmpmm'))  );
	

	#my ($otherformtype) = do 
	#{	
	#	if ($formtype eq 'simple') {'advanced'}
	#	elsif ($formtype eq 'advanced') {'simple'}
	#	else {''};
	#};

	my (%modes);
	%modes=(
					'one_strand'=>'generic motifs on one strand',
					'two_strands'=>'generic motifs on both strands',
					'palindromes'=>'palindromes',
					'one_strand_rep'=>'repeats on one strand',
					'two_strands_rep'=>'repeats on both strands',
					);
	my ($arga_txt);

	if ($is_example)
	{
$arga_txt=">argR\natgtttctcaataacgaaatttgataaaatcccgctctttcataacattatttcagccttcttcagggctgactgtttgcataaaaattcatctgtatgcacaataatgttgtatcaaccaccatatcgggtgacttatgcgaagctcggctaagcaagaagaactagttaaagcatttaaagcattacttaaagaagag\n".
">argA\nttcacggcattactgataaaaaagtcgctctcgcataaaatttacacttgcaccctgcgaaaaaacagaataaaaatacactaatttcgaataatcatgcaaagaggtgtgccgtggtaaaggaacgtaaaaccgagttggtcgagggattccgccattcggttccctatatcaatacccaccggggaaaaacgtttgtc\n".
">argCBH\nataaatggcggtaatttgtttttcattgttgacacacctctggtcatgatagtatcaatattcatgcagtatttatgaataaaaatacactaacgttgagcgtaataaaacccaccagccgtaaggtgaatgttttacgtttaacctggcaaccagacataagaaggtgaatagccccgatgttgaatacgctgattgtg\n".
">argD\ngtttgttcaattgccatctcatgatcaccctgttacgcataaacaaatgtgaaattataaccacaaaatatgcataaaaaatcactaaatggcaatcagaaatcagcgatgcaggaaattagccagcagttgatgtccttgttcgctaagaatactttctggatggaactgcacaccttccagatcccactggcgatggcgaatccccataatctctcgggtttcgctccaggccgtcac\n".
">argE\nataaatggcggtaatttgtttttcattgttgacacacctctggtcatgatagtatcaatattcatgcagtatttatgaataaaaatacactaacgttgagcgtaataaaacccaccagccgtaaggtgaatgttttacgtttaacctggcaaccagacataagaaggtgaatagccccgatgttgaatacgctgattgtg\n".
">argF\nagaaagtgttttttgtataaatcggacattttatcctcgcatggcgaacgccacttattgaattaaaattcactttatatgtgtaattattcatttgcaaccccatttcacaattctttcttacaaaggtggaggcaaacccgtccgtgtgtgaaaataatcgtatctgcctccgattctctgcagaagcagaaagacat\n".
">argG\naacgtttattgctaatcatgtgaatgaatatccagttcactttcatttgttgaatacttttgccttctcctgctctcccttaagcgcattattttacaaaaaacacactaaactcttcctgtctccgataaaagatgattaaatgaaaactcatttattttgcataaaaattcagtgagagcggaaatccaggctcatcacgacgagctaccaggctgctccaccccgcgcctgaaacgtggcaaattctactcgttttgggtaaaaaatgcaaatactgctgggatttggtgtaccgagacgggacgtaaaatctgcaggcattatagtgatccacgccacattttgtcaacgtttattgctaatcatgtgaatgaatatccagttcactttcatttgt\n".
">argI\natgcttatgataaaacccggacatagatccctcctgtggctaacgcctcaatgaattaaaattcaatttatatggatgattattcatttgcaagtctaaagcataaatctttgtcacaaaggtggaggcaatgtcagtggtgtgtgacaataagagtatcggcaggacattaagaggaatgagccatggcaaacccggaa\n".
">carAB\nttgtgaattaatatgcaaataaagtgagtgaatattctctggagggtgttttgattaagtcagcgctattggttctggaagacggaacccagtttcacggtcgggccataggggcaacaggttcggcggttggggaagtcgttttcaatacttcaatgaccggttatcaagaaatcctcactgatccttcctattctcgt";
	}

	my ($takemaskedfasta,$prev_id)=($form->param('takemaskedfasta'),$form->param('prev_id'));

	my ($mfasta);
	
	if ((defined $takemaskedfasta) && $takemaskedfasta eq 'yes')
	{
		my $dirname=$WORKDIR.'/gibbslfm.'.$prev_id;
		my $fname="";
		open IN, "$dirname/result";
		while(<IN>)
		{
			if (/\$(.*\S)\s*--/ || /\$(.*\S)\s*$/)
			{
				$fname=$1;
				last;
			}
			if (/\$\s*--/ || /\$\s*$/)
			{
				$fname="";
				last;
			}
		}
		close IN;

		if (-s $WORKDIR.'/gibbslfm.'.$prev_id.'/FastA.out.masked')
		{
			open ( MFASTA, $WORKDIR.'/gibbslfm.'.$prev_id.'/FastA.out.masked');
			while (<MFASTA>)
			{
				$mfasta=$mfasta.$_;
			}
			close MFASTA;
			$form->param(-name=>'fasta_text',-value=>$mfasta);
			if ($fname=~/(.+)\.masked-(\d*)\.(\w+)/)
			{
				$fname=$1.'.masked-'.($2+1).$3;
			}
			else
			{
				$fname=~/(.+)\.(\w+)/;
				$fname=$1.'.masked-1.'.$2;
			}
			$form->param(-name=>'fasta_file_name',-value=>$fname);	
		}
		else
		{
			print	$form->header, 
						$form->start_html(-title=>'Wait for task '.$prev_id.' completion.',
						-bgcolor=>'white'),
						$form->h3('Please, wait for task with id '.$prev_id.' (Fasta file name '.$fname.').'),
						$form->end_html;									
			return;
		}
	}
	else
	{$takemaskedfasta='no'}
#takemaskedfasta=yes&prev_id=id
	print $form->header,
				$form->start_html(-title=>ucfirst($formtype).' form for SeSiMCMC data&parameters input.',-script=>{-language=>'JavaScript',-code=>$Script},
				-bgcolor=>'white');
				&print_header_links;

	print	$form->start_multipart_form(-action=>$url,-method=>'post',-name=>'MainInputForm'),
				$form->h2(ucfirst($formtype).' form for SeSiMCMC data&amp;parameters input'),
				$form->a({href=>$url.'?action=cat'.
												"&source=readme&format80=yes"},
										"I want to read the readme file explaining the options." ),
				$form->hr,"\n",
				$form->b('I want to recall results for id: '),
				$form->textfield(	-name=>'id',-size=>8),
				$form->submit(-name=>'action',-value=>'<< Recall'),
				$form->hr,"\n",
				'<b><input type="submit" style="height: 50px" name="action", value="Run sampler &gt;&gt;"</b>',
				#$form->b($form->submit(-name=>'action',-value=>'Run sampler >>')),
				"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
				#$form->b($form->submit(-name=>'action',-value=>'Completion example')),
				'<b><input type="submit" style="height: 50px" name="action", value="Completion example"</b>',
				do {
					if ($takemaskedfasta ne 'yes')
					{
						$form->h3('Type the name or browse your FastA file here'). 
						$form->filefield(-name=>'fasta_file',-size=>50,-onChange=>'form.elements[\'fasta_file_name\'].value=form.elements[\'fasta_file\'].value').$form->br;}
				},
				$form->h3('The FastA filename to be referred in output:'),
				do {
					if ($is_example) {
						$form->textfield(-name=>'fasta_file_name',-size=>50,-value=>'ArgR.dnc',-override=>1);
					} 
					else { 
						$form->textfield(-name=>'fasta_file_name',-size=>50);
					} 
				},
				$form->h3('Type or paste your FastA data here'),
				do {
					if ($is_example) {
						$form->textarea(-name=>'fasta_text',
														-rows=>10,
														-columns=>50,
														-value=>"$arga_txt",
														-override=>1),
					} 
					else { 
						$form->textarea(-name=>'fasta_text',
														-rows=>10,
														-columns=>50),
					} 
				},
				$form->br,
#            1234567890123456789012345678901234567890
				'<table width=600>',
				'<tbody>',
				'<tr></tr>',
				'<tr><td>',
				'<TT>We are looking for:</TT>',
				'</td><td colspan=3>',
			 	$form->popup_menu(
							 -name=>'mode',
							 -values=>[$modes{'one_strand'},
							 					 $modes{'two_strands'},
												 $modes{'palindromes'},
												 $modes{'one_strand_rep'},
												 $modes{'two_strands_rep'}],
							 -default=>$modes{'two_strands'}),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
				'<TT>Start motif length:</TT>', 
				'</td><td>',
				$form->textfield(	-name=>'motif_length',-size=>8,
													-value=>'default',
										    	-onFocus=>'SelectIfDefault(name)',
													-onBlur=>'LeaveIntField(name)'),
				'</td><td></td><td align=right>',
				'<TT>Adjust motif length:</TT>',
				'</td><td>',
				$form->checkbox(-name=>'adjust_motif_length',
                        -checked=>'checked',
                        -value=>'yes',
												-label=>''),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
				'<TT>Minimal motif length:</TT>', 
				'</td><td>',
				$form->textfield(	-name=>'minimal_motif_length',-size=>8,
													-value=>'default',
										    	-onFocus=>'SelectIfDefault(name)',
													-onBlur=>'LeaveIntField(name)'
											),
				'</td><td>',
				'</td><td align=right>',
				'<TT>Slow optimisation:</TT>',
				'</td><td>',
				$form->checkbox(-name=>'slow_optimisation',
                        -value=>'yes',
												-label=>''),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
				'<TT>Maximal motif length:</TT>', 
				'</td><td>',
				$form->textfield(	-name=>'maximal_motif_length',
									-size=>8,
									-value=>'default',
									-onFocus=>'SelectIfDefault(name)',
									-onBlur=>'LeaveIntField(name)'
											),
				'</td><td>',
				'</td><td align=right>',
				'<TT>Possibly spaced motif:</TT>',
				'</td><td>',
				$form->checkbox(-name=>'spaced_motif',
                        -value=>'yes',
						-label=>''),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
				'<TT>Motif absence prior (probability):</TT>',
				'</td><td>';
			    
				#here, cgi.pm textfield() function became crasy. Was workaround.	
				print '<input type="text" name="motif_absence_prior" value="',
				do {
					if($formtype eq 'dmmpmm')
					{"0";}
					else
					{"default";}	
				},
				'" size="8" onfocus="SelectIfDefault(name)" onblur="LeaveProbabilityField(name)" /></td><td>',
				'</td><td align=right>',
				'<TT>Output FastA with masked sites:</TT>',
				'</td><td>',
				$form->checkbox(-name=>'masked_output',
                        -value=>'yes',
												-label=>''),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td colspan=4>',
				'<TT>Retrieve additional sites that are not worse than the profile sites',
				'</td><td>',
				$form->checkbox(-name=>'retrieve_other_sites',
                        						-value=>'yes',
												-checked=>'checked',
												-label=>''),
				'</td></tr>'; 
	

		print 	'<tr></tr>',
		      	'<tr><td>',
				'<TT>Background model is common (static):</TT>',
				'</td><td>';
				if ($formtype eq 'dmmpmm')
				{
					print '<input type="checkbox" name="common_background" value="yes" checked="checked"/>',
					'</td><td></td><td>',
				'<TT>Add poly-n 100 base flanks to each sequence head and tail:</TT>',
				'</td><td><input type="checkbox" name="poly_n_flanks" value="yes" checked="checked"/></td></tr>';
				}
				else
				{
					print '<input type="checkbox" name="common_background" value="yes" />',
					'</td><td></td></tr>';
				};

	if ($formtype eq 'dmmpmm')
	{
		print	'<tr><td colspan=3>',
				'<TT>Caps in input sequences are anchors for the motif sites. It is useful not to get poly-n in sites. Do not forget to capitalise some letters in source sequences.</TT>',
				'</td><td>',
				$form->checkbox(-name=>'obligatory_caps',
                	       -value=>'yes',
                   			-label=>''),
				'</td></tr>',
			
	};
	
	
	if ($formtype eq 'advanced')
	{
		print	'<tr></tr><td>',
				'<TT>Caps in input sequences are anchors for the motif sites:</TT>',
				'</td><td>',
				$form->checkbox(-name=>'obligatory_caps',
                	       -value=>'yes',
							-label=>''),
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
		       	'<TT>Random seed #1:           </TT>',	
				'</td><td>',
			    $form->textfield (-name=>'random_seed_1',-size=>8,
														-value=>'default',
										    		-onFocus=>'SelectIfDefault(name)',
														-onBlur=>'LeaveIntField(name)'),
 				'</td><td>',
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
		       	'<TT>Random seed #2:           </TT>',	
				'</td><td>',
			    $form->textfield (-name=>'random_seed_2',-size=>8,
														-value=>'default',
										    		-onFocus=>'SelectIfDefault(name)',
														-onBlur=>'LeaveIntField(name)'),
				'</td><td>',
				'</td></tr>',
				'<tr></tr>',
				'<tr><td>',
		       	'<TT>Chains to try:           </TT>',	
				'</td><td>',
			    $form->textfield (-name=>'chains_to_try',-size=>8,
														-value=>'default',
										    		-onFocus=>'SelectIfDefault(name)',
														-onBlur=>'LeaveIntField(name)'),
				'</td><td>',
				'</td></tr>',
				'<tr></tr>',
#				'<tr><td nobreak>',
#				'<TT>"Annealed" test is:</TT>',
#					'</td><td>',
#			 	  $form->popup_menu(
#							 -name=>'annealing',
#							 -values=>['strong','usual','weak']),
#				  '</td></tr>',
					 '<tr><td nobreak>',
				   '<TT>A maximum holds for</TT>',
					 '</td><td>',
			 	   $form->popup_menu(
							 -name=>'annealings_number_maximum_is_to_be_global_for',
							 -values=>['5','10','20','50']),
					 '</td><td>',
					 '</td><td>',
				   '<TT>annealing times to be global.</TT>',
				   '</td></tr>',
				   '<tr></tr>',
					 '<tr><td>',
					 '<TT>Part of every found site to mask on output:</TT>',
					 '</td><td>',
					 $form->textfield(	-name=>'masked_part',-size=>8,
														-value=>'default',
														-onFocus=>'SelectIfDefault(name)',
														-onBlur=>'LeaveProbabilityField(name)'
												),
					 '</td></tr>',
					 '</table>',
					 '<table width=400>',
				   '<tr></tr>',
					 '<tr><td colspan=4>',
					 '<TT>Set the background manually:</TT>&nbsp',
					 $form->checkbox(-name=>'manual_background',
													-value=>'yes',
													-label=>'',
													-onClick=>'SyncBackgrGrey()'),
					 '</td></tr>',
				   '<tr></tr>',
					 '<tr><td>',
					 '<TT>A:</TT>&nbsp',
  		     $form->textfield (-name=>'background_A',-size=>8,
														-value=>'auto',
														-onBlur=>'LeaveBackgroundField(name)'),
					 '</td><td>',
					 '<TT>T:</TT>&nbsp',
			     $form->textfield (-name=>'background_T',-size=>8,
														-value=>'auto',
														-onBlur=>'LeaveBackgroundField(name)'),
					 '</td><td>',
					 '<TT>G:</TT>&nbsp',
			     $form->textfield (-name=>'background_G',-size=>8,
														-value=>'auto',
														-onBlur=>'LeaveBackgroundField(name)'),
					 '</td><td>',
					 '<TT>C:</TT>&nbsp',
			     $form->textfield (-name=>'background_C',-size=>8,
														-value=>'auto',
														-onBlur=>'LeaveBackgroundField(name)'),
					 '</td></tr>',
	};
	print '</table>',"\n";

	if ($formtype eq 'simple' or $formtype eq "DMMPMM")
	{
		print	$form->hidden (-name=>'chains_to_try',
												  -value=>'default');
#		print	$form->hidden (-name=>'annealing',
#												  -value=>'usual');
		print	$form->hidden (-name=>'annealings_number_maximum_is_to_be_global_for',
												  -value=>'10');
		print	$form->hidden (-name=>'random_seed_1', -value=>'default');
		print	$form->hidden (-name=>'random_seed_2', -value=>'default');
		print	$form->hidden (-name=>'masked_part', -value=>'default');
	}
#we use alias "annealing=usual for 
#cycles_with_minor_change_to_say_cold=4
#minor_change_level=0.5
#we use alias "annealing=weak for 
#cycles_with_minor_change_to_say_cold=2
#minor_change_level=0.1
#we use alias "annealing=strong for 
#cycles_with_minor_change_to_say_cold=5
#minor_change_level=0.75

#and we alias modes as described in %modes
# all other, we use tags from config.
#for boolean tags, we store "yes" if on
#and do not store it or set "no" if no. 
	
	#local_step_cycles_between_adjustments=10
	print	$form->hidden (-name=>'local_step_cycles_between_adjustments', -value=>'default'),
  #adjustments_during_annealing=1
				$form->hidden (-name=>'adjustments_during_annealing',
											 -value=>'yes'),
	#steps_number_maximum_is_to_be_global_for=0
				$form->hidden (-name=>'steps_number_maximum_is_to_be_global_for',
									-value=>'default'),
	#chain_fails_after=3
				$form->hidden (-name=>'chain_fails_after',
												-value=>'default'),"\n";
	print $form->hr;
	print '<table border="0" columns="3"><tr>';
	if ($formtype eq 'simple')
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Advanced \\/'),'</td>';
	}
	elsif ($formtype eq 'advanced')
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Simple /\\'),'</td>';
	}
	else #dmmpmm
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Simple /\\'),'</td>';
	}
	
#	print '<td></td><td rowspan="3">',$form->b($form->submit(-name=>'action', -height=>"200" ,-value=>'Run sampler >>')),'</td>';
	print '<td></td><td rowspan="3"><b><input type="submit" style="height: 50px" name="action", value="Run sampler &gt;&gt;"</b></td>';
	
	print '</tr><tr></tr><tr>';		

	if ($formtype eq 'simple')
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Footprint search as in DMMPMM '),'</td>';
	}
	elsif ($formtype eq 'advanced')
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Footprint search as in DMMPMM '),'</td>';
	}
	else #dmmpmm
	{
		print	'<td>',$form->submit(-name=>'action',-value=>'Advanced \\/'),'</td>';
	}
	print			'</tr></table>';
				
	print		$form->endform;
	my $goole_anal_script=
	'<script type="text/javascript">'."\n".
	'var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");'."\n".
	'document.write(unescape("%3Cscript src='."\'".'" + gaJsHost + "google-analytics.com/ga.js'."\'".
	' type='."\'".'text/javascript'."\'".'%3E%3C/script%3E"));'."\n".
	'</script>'."\n".
	'<script type="text/javascript">'."\n".
	'try {'."\n".
	'var pageTracker = _gat._getTracker("UA-7441172-1");'."\n".
	'pageTracker._trackPageview();'."\n".
	'} catch(err) {}</script>'."\n";
	print $goole_anal_script;
	print	$form->end_html;
}

sub run_sampler
{
	my ($cgi) = @_;
	my ($url)=$cgi->url();

	$url=~s/\?.*//;
	my ($id,$dirname,$fasta_deliver,$fasta_text,$key,$value,$runerror,$fasta_file_name);

	my $cmdstr = "mktemp -d $WORKDIR/gibbslfm.XXXXXX";
	$dirname=`$cmdstr`;
	my $chmodstr= "chmod 775 $dirname";
	`$chmodstr`;
	chomp $dirname;
	$id = substr $dirname, -6;
	$cgi->param(-name=>'id',-value=>$id);

	$fasta_text=$cgi->param('fasta_text');
	$fasta_deliver='textarea' if (defined $fasta_text && $fasta_text);

	my $filename=$cgi->param('fasta_file');
	$fasta_deliver='upload' if (defined $filename && $filename);
	
	$fasta_file_name=$cgi->param('fasta_file_name');
	$fasta_file_name=$filename 
		if (!defined $fasta_file_name && $fasta_deliver eq 'fasta_file');
	
	open (FASTA,">$dirname/FastA.seq") or $fasta_deliver='error';

	if ($fasta_deliver eq 'upload')
	{
		if (!$filename)
		{
			print $cgi->header,
						$cgi->start_html(-title=>'An error occured on target file opening.',
						-bgcolor=>'white'),
						$cgi->h2("Empty file name."),
						$cgi->end_html;
						return;
		}
		while (<$filename>)
		{
			print FASTA;
		}
		close FASTA;
		$fasta_file_name=$filename;
	}
	elsif ($fasta_deliver eq 'textarea')
	{
	  print FASTA $cgi->param('fasta_text');
		close FASTA;
	}

	elsif ($fasta_deliver eq 'error')
	{
		print $cgi->header,
					$cgi->start_html(-title=>'An error occured on target file opening.',
					-bgcolor=>'white'),
					$cgi->h2("I cannot open file $dirname/FastA.seq."),
  				$cgi->end_html;
					return;
	}
	else
	{
		print $cgi->header,
					$cgi->start_html(-title=>'An error occured in the script.',
					-bgcolor=>'white'),
					$cgi->h2('Something is wrong. I don\'t know what to do without FastA data.'),
  				$cgi->end_html;
					close FASTA;
					system "rm -rf $dirname";
					return;
	}
	
	my (%modes_code);
	%modes_code=(
					'generic motifs on one strand'=>'one_strand',
					'generic motifs on both strands'=>'two_strands',
					'palindromes'=>'palindromes',
					'repeats on one strand'=>'one_strand_rep',
					'repeats on both strands'=>'two_strands_rep',
					);
					
	my ($masked_output_key)=('');
	if ($cgi->param('masked_output') eq 'yes') 
	{
		$masked_output_key="--output-fasta $dirname/FastA.out.masked ";
		$masked_output_key=$masked_output_key.'--masked-part '.
				$cgi->param('masked_part').' ' if $cgi->param('masked_part') ne 'default';
	}

#it is for future start...
	my ($run_command) = "nice $sampler_name --cgi --html --output-tables --ip $ENV{REMOTE_ADDR} --log $WORKDIR/logfile --input-name \'$fasta_file_name\' --id $id -i $dirname/FastA.seq --config-file $dirname/config --output-file $dirname/result  $masked_output_key>$dirname/error 2>&1";

#it is to generate default configuration
	my ($conf_command) = "$sampler_name -q -i $dirname/FastA.seq -w --config-file $dirname/config >$dirname/config.out 2>$dirname/error && mv $dirname/config.out $dirname/config";

#results file plug preparation...
	open(RESULT,">$dirname/result");
	print RESULT $cgi->start_html(-title=>'"Waiting for results..."',
													-head=>$cgi->meta(
														{ 
															-http_equiv=>'refresh',
															-content=>'5'
														}
													),-bgcolor=>'white'),
				'<!-- $'.$fasta_file_name.'-->'."\n",
				$cgi->h3('The task with id '.$id.' ('.$fasta_file_name.') is working...'),"\n",
				'<!-- '.$run_command.' -->'."\n",
				$cgi->end_html;
	close(RESULT);
#the config file preparation..
	open (CONF,">$dirname/config");
	my ($modestri);
	$modestri=$modes_code{$cgi->param('mode')};
	print CONF 'mode=', 
		do {
			if ( $modestri=~/one_strand/ ) {"one_strand";}
			else {"two_strands";}
			},"\n";
	print CONF 'symmetry=', 
		do {
			if ( $modestri	eq "palindromes" ) 
				{"palindromes";}
			else 
			{ 
				if ( $modestri=~/rep/ ) {"repeats";}
				else
				{"no";}
			}
		},"\n";
	print CONF 'motif_absence_prior=',$cgi->param('motif_absence_prior'),"\n" 
		if ( defined $cgi->param('motif_absence_prior') && 
					$cgi->param('motif_absence_prior') ne 'default');
	print CONF 'adjust_motif_length=',
						 	do {
								if ($cgi->param('adjust_motif_length') eq 'yes') 
									{'yes';}
								else {'no'};
							},"\n",
						 'spaced_motif=',
						 	do {
								if ($cgi->param('spaced_motif') eq 'yes') 
									{'yes';}
								else {'no'};
							},"\n";
							
	print CONF 'motif_length=',$cgi->param('motif_length'),"\n" 
		if ( defined $cgi->param('motif_length') && 
					$cgi->param('motif_length') ne 'default');
					
	print CONF 'maximal_motif_length=',$cgi->param('maximal_motif_length'),"\n" 
		if ( defined $cgi->param('maximal_motif_length') && 
					$cgi->param('maximal_motif_length') ne 'default');
	print CONF 'minimal_motif_length=',$cgi->param('minimal_motif_length'),"\n" 
		if ( defined $cgi->param('minimal_motif_length') && 
					$cgi->param('minimal_motif_length') ne 'default');
					
					
	print CONF 'retrieve_other_sites=',
						 	do {
								if ($cgi->param('retrieve_other_sites') eq 'yes') 
									{'smart';}
								else {'no'};
							},"\n";
							
	#print CONF 'retrieve_other_sites_threshold=',$cgi->param('retrieve_other_sites_threshold'),"\n"
	#	if (defined $cgi->param('retrieve_other_sites_threshold'));

	print CONF 'local_step_cycles_between_adjustments=',
			$cgi->param('local_step_cycles_between_adjustments'),"\n" 
		if ( defined $cgi->param('local_step_cycles_between_adjustments') && 
					$cgi->param('local_step_cycles_between_adjustments') ne 'default');

	print CONF 'adjustments_during_annealing=',
							do {
								if ($cgi->param('adjustments_during_annealing') eq 'yes') 
									{'yes';}
								else {'no'};
							#},
							#"\n",
						 	#do {
						 	#	if($cgi->param('annealing') eq 'usual')
							#	{
							#		"cycles_with_minor_change_to_say_cold=5\n".
							#		"minor_change_level=0.7";
							#	}
							#	elsif ($cgi->param('annealing') eq 'weak') 
							#	{
							#		"cycles_with_minor_change_to_say_cold=2\n".
							#		"minor_change_level=0.3";
							#	}
							#	elsif ($cgi->param('annealing') eq 'strong') 
							#	{
							#		"cycles_with_minor_change_to_say_cold=5\n".
							#		"minor_change_level=0.75";
							#	}
	            },"\n";

	print CONF 'annealings_number_maximum_is_to_be_global_for=',
			$cgi->param('annealings_number_maximum_is_to_be_global_for'),"\n" 
				if ( defined 
					$cgi->param('annealings_number_maximum_is_to_be_global_for') && 
					$cgi->param('annealings_number_maximum_is_to_be_global_for') ne
									'default');

	print CONF 'steps_number_maximum_is_to_be_global_for=',
			$cgi->param('steps_number_maximum_is_to_be_global_for'),"\n" 
				if ( defined $cgi->param('steps_number_maximum_is_to_be_global_for') && 
					$cgi->param('steps_number_maximum_is_to_be_global_for') ne 'default');

	print CONF 'slow_optimisation=',
						 	do {
								if ($cgi->param('slow_optimisation') eq 'yes') 
									{'yes'}
								else {'no'}
							},"\n";
							
	print CONF 'chain_fails_after=',$cgi->param('chain_fails_after'),"\n" 
		if ( defined $cgi->param('chain_fails_after') && 
					$cgi->param('chain_fails_after') ne 'default');
					
	print CONF 'chains_to_try=',$cgi->param('chains_to_try'),"\n" 
		if ( defined $cgi->param('chains_to_try') && 
					$cgi->param('chains_to_try') ne 'default');
					
	print CONF 'input_caps=',
						 	do {
								if ($cgi->param('obligatory_caps') eq 'yes') 
									{'obligatory'}
								else {'neglect'}
							},"\n";

	print CONF 'common_background=',
						 	do {
								if ($cgi->param('common_background') eq 'yes') 
									{'yes'}
								else {'no'}
							},"\n";

	print CONF 'random_seed_1=',$cgi->param('random_seed_1'),"\n" 
		if ( defined $cgi->param('random_seed_1') && 
					$cgi->param('random_seed_1') ne 'default');
					
	print CONF 'random_seed_2=',$cgi->param('random_seed_2'),"\n" 
		if ( defined $cgi->param('random_seed_2') && 
					$cgi->param('random_seed_2') ne 'default');

	if ($cgi->param('manual_background') eq 'yes')
	{
		print CONF 'background_A=',$cgi->param('background_A'),"\n"; 
		print CONF 'background_T=',$cgi->param('background_T'),"\n"; 
		print CONF 'background_G=',$cgi->param('background_G'),"\n"; 
		print CONF 'background_C=',$cgi->param('background_C'),"\n"; 
	}
					
  close CONF;
  
	#preparing flanked FASTA
	if ($cgi->param('poly_n_flanks') eq 'yes')
	{
		my $poly_n="";
		for (my $pos=1; $pos<=100; $pos++)
		{
			$poly_n='n'.$poly_n;
		}
		my ($seqname,$sequence);
		open (FASTA_IN,"$dirname/FastA.seq");
		open (FASTA_OUT,">$dirname/FastA.flanked.seq");
		while(<FASTA_IN>)
		{
			chomp;
			s/\r//g;
			 if (/^\s*#/) {print FASTA_OUT ; next;}; # omits empty line
			if (!defined $seqname)
			{
			#if the name is defined then we wait for it
				if (/^\s*$/) {print FASTA_OUT ; next;} ; # omits empty line
				if ((substr $_,0,1) eq ">") # a name; a new line!
				{
					$seqname=substr $_,1;
					$sequence="";
				}
			}
			else  #we are reading something now
			{
				if (/^\s*$/) #empty line
				{
					print FASTA_OUT '>',$seqname,"\n";
					print FASTA_OUT $poly_n.$sequence.$poly_n,"\n";
					$sequence=undef;
					$seqname=undef;
				}
				elsif ((substr $_,0,1) eq ">") # a name; a new line!
				{
					print FASTA_OUT '>',$seqname,"\n";
					print FASTA_OUT $poly_n.$sequence.$poly_n,"\n";
					$seqname=substr $_,1;
					$sequence="";
				}
				else #elongation of sequence
				{
					$sequence=$sequence.$_;
				}
			}
		}
		if (defined $seqname) #do not loose the last one!
		{
			print FASTA_OUT '>',$seqname,"\n";
			print FASTA_OUT $poly_n.$sequence.$poly_n,"\n";
			$sequence=undef;
			$seqname=undef;
		}
		close FASTA_OUT;
		close FASTA_IN;
	}
 	
#making the page redirection for index; 

	print $page->header,
  			$page->start_html(-title=>'"Creating index..."', 
													-head=>$page->meta(
														{ 
															-http_equiv=>'refresh',
															-content=>
																	'0; url='.$url.'?action=index'.'&id='.$id
														}
													)),
				$page->end_html;


	$runerror=system($conf_command);
#preparing a command and start the sampler
	if (!$runerror)
	{
		# do smth like: $runerror=system("$run_command");
		my $kid = fork();
		unless (defined $kid) {
			$runerror = "Fatal error. Seems we can't fork";
		} else {
			if ($kid) { # parent
				close STDOUT;
				my ($ret);
	#			warn "Wait start..";
				do {
					$kid = waitpid(-1,0);
				} until $kid == -1;
	#			$runerror=$? if $?;
	#			warn "Wait returned $ret";
			} else { # child
				close STDOUT;
				exec($run_command) or $runerror = "exec: $?";
				warn "Wow! We're still alive...";
			}
		}
	}
	if ($runerror)
	{
		open(RESULT,">$dirname/result");
		print RESULT $cgi->start_html(-title=>'"Something is wrong..."',
		             -bgcolor=>'white'), 
					'<!-- $'.$fasta_file_name.'-->'."\n",
					$cgi->h3('The task with id '.$id.' ('.$fasta_file_name.') was finished with an error '. $runerror),
					$cgi->a({href=>$url.'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/error'},'Try to look at diagnostics'),"\n<br>",
					$cgi->end_html;
		close(RESULT);
	}
	#else
	#{
	#	warn "ren";
  #	rename ("$dirname/run_results","$dirname/result");
	#}
	
}

sub cat
{
	my ($outpage) = @_;
	my ($fname,$ifformat80,$html);
	my $filename = $outpage->param('source');
	my $allowed=0;
	$allowed=1 if (($filename eq 'readme') or ($filename=~/^$WORKDIR/));
	$ifformat80 = defined $outpage->param('format80');
	open IN, $filename;
	my $header_printed;
	if (! $allowed) #non-allowed file
	{
		print $outpage->header,
		      $outpage->start_html(-title=>$filename,-bgcolor=>'white'),
					">>You don't have permissions to cat the file.<<",
					$outpage->end_html;
		return;
	}
	while (<IN>){
		if (! defined $html)
		{
			if (!$header_printed)
			{
				print $outpage->header;
				$header_printed=1;
			}
			if (/\S/) #a nonempty line
			{
				if (/^\s*<\s*\?\s*xml/i)
				{
					print; #we do not eat xml tag, still we do not decide what to do
				}
				elsif (/^\s*<\s*\!\s*DOCTYPE\s*HTML/i)
				{
					$html=1;
				}
				elsif (/<\s*HTML/i)
				{
					$html=1;
				}
				else
				{
					$html=0;
					print $outpage->start_html(-title=>$filename,-bgcolor=>'white'),
								"<pre width=80>\n";
				}
			}
			else
			{
				next; 	
			}
		}
		if (defined $html) # it is not "else", possibly it is defined just now
		{
			if ($html)
			{
				print;
			}
			else
			{
				s/</&lt;/g;
				s/>/&gt;/g;
				if ($ifformat80)
				{
					my $position=0;
					my @line=split(/ /);
					foreach my $token (@line)
					{
						$position+=(1 + length($token));
						if ($position>80)
						{
							$position=1 + length($token);
							print "\n";
						}
						print " ",$token;
					}
				}
				else
				{
					print;
				}
			}
		}
	}
	close IN;
	
	if (defined $html and !$html)
	{
		print "</pre>", $outpage->end_html;
	};
	
	if (! defined $html) #empty file
	{
		print $outpage->header,
		      $outpage->start_html(-title=>$filename,-bgcolor=>'white'),
					">>The file is empty<<",
					$outpage->end_html;
	}
}

sub makerindexredir
{
	my ($page) = @_;
	my ($url)=$page->url();
	$url=~s/\?.*//;
	my ($id)=$page->param('id');
	print $page->header,
  			$page->start_html(-title=>'"Creating index..."', 
													-head=>$page->meta(
														{ 
															-http_equiv=>'refresh',
															-content=>
																	'0; url='.$url.'?action=rindex'.'&id='.$id
														}
													)),
				$page->end_html;
	
}

sub makeindex
{
	my ($outpage,$isrecall)= @_;
	my ($id,$dirname,$fasta_file_name);
	$id=$outpage->param('id');
	$dirname=$WORKDIR.'/gibbslfm.'.$id;
	if (! -d $dirname)
	{
		print	$outpage->header, 
					$outpage->start_html(-title=>'The SeSiMCMC task '.$id.' information.',
					-bgcolor=>'white');
		&print_header_links;
		print	$outpage->h3('There is no information about the task with id "'.$id.'", sorry.'),
					$outpage->hr,
					$outpage->a({href=>$outpage->url().'?action=form'},'Go back to data input page'),
					"\n<br>",
					$outpage->end_html;									
		return;
	}
	open IN, "$dirname/result";
	while(<IN>)
	{
		if (/\$(.*\S)\s*--/ || /\$(.*\S)\s*$/)
		{
			$fasta_file_name=$1;
			last;
		}
		if (/\$\s*--/ || /\$\s*$/)
		{
			$fasta_file_name="";
			last;
		}
	}
	close IN;
	print	$outpage->header;
	&print_header_links;
	print	$outpage->start_html(-title=>'The SeSiMCMC task '.$id.' information.',
				-bgcolor=>'white');
				if (defined $isrecall && ! $isrecall)
				{
					print $outpage->h3('Your task identificator is '.$id.', remember it if you want to access the data later.');
				}
				else
				{
					print $outpage->h3('The task with identificator: '.$id.' is recalled.');
				};
				print $outpage->h3('The FastA file name is "'.$fasta_file_name.'".'),
				$outpage->hr,
				$outpage->h3($outpage->a({href=>$outpage->url().'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/result'},'Results')),"\n",
				$outpage->hr,
				$outpage->a({href=>$outpage->url().'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/FastA.seq'},'FastA data'),"\n<br>",
				$outpage->a({href=>$outpage->url().'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/config'},'Parameters (the configuration file)'),"\n<br>",
				$outpage->a({href=>$outpage->url().'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/error'},'Diagnostics'),"\n<br>",
				do {
					if ( -f $WORKDIR.'/gibbslfm.'.$id.'/FastA.out.masked') 
						{
							$outpage->a({href=>$outpage->url().'?action=cat&source='.$WORKDIR.'/gibbslfm.'.$id.'/FastA.out.masked'},'FastA with masked sites')."&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp".$outpage->a({href=>$outpage->url().'?action=form&takemaskedfasta=yes&prev_id='.$id},'Start the sampler with the masked data')."\n<br>";
						}
				},
				$outpage->hr,
				$outpage->a({href=>$outpage->url().'?action=form'},'Go back to data input page'),"\n<br>";
								
				
	print $outpage->end_html;									
}

sub print_version
{
	my ($outpage) = @_;
	print	$outpage->header, 
				$outpage->start_html('SeSiMCMC (Sequence Similarities MCMC) version information',
				-bgcolor=>'white'),
	      "The version of current SeSiMCMC executable follows:<BR>\n",
	      '<HR>',
				'<TT>';
	open (INFILE,"$sampler_name --version |");
	while (<INFILE>){
		print; 
		print '<BR>'; 
		print "\n"; 
	}
  print '</TT>';
	print '<HR>';
	print $outpage->end_html;
} 

sub print_header_links
{
print '<p><a href=http://bioinform.genetika.ru>Bioinformatics Laboratory, State Scientific Centre GOSNIIGenetika, Moscow</a></p>';
print '<p><a href=',$SeSiMCMC_home,'>The SeSiMCMC home page</a></p><hr>';
}
