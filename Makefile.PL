use ExtUtils::MakeMaker;

WriteMakefile(
	NAME              => 'Stats::LikeR',
	VERSION_FROM      => 'lib/Stats/LikeR.pm', 
	PREREQ_PM         => {},
	CCFLAGS           => '-std=c99',
	LIBS              => ['-lm'], 
	DEFINE            => '', 
	INC               => '-I.', 
    # Primary optimization control
	OPTIMIZE          => '-O4 -Wunused',
    # For linking step as well
);
