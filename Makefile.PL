use ExtUtils::MakeMaker;

WriteMakefile(
	NAME              => 'Stats::LikeR',
	VERSION_FROM      => 'lib/Stats/LikeR.pm', 
	PREREQ_PM         => {}, 
	LIBS              => ['-lm'], 
	DEFINE            => '', 
	INC               => '-I.', 
);
