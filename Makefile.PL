use ExtUtils::MakeMaker;

WriteMakefile(
	NAME              => 'stats',
	VERSION_FROM      => 'lib/stats.pm', 
	PREREQ_PM         => {}, 
	LIBS              => ['-lm'], 
	DEFINE            => '', 
	INC               => '-I.', 
);
