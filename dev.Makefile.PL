use ExtUtils::MakeMaker;

WriteMakefile(
	NAME              => 'Stats::LikeR',
	VERSION_FROM      => 'lib/Stats/LikeR.pm', 
	PREREQ_PM         => {},
	CCFLAGS           => '-std=gnu99',
	LIBS              => ['-lm'], 
	DEFINE            => '', 
	INC               => '-I.', 
);
