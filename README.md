# IPCC2024Repo
This a repository aims to reproduce the results shown in a paper submitted to ICPP 2024

Three artifacts are used in this paper. Based on the source-to-source transformation compiler ROSE, we developed the toolkit “UQIntrinsicTransform” as an implementation of the proposed uncertainty propagation detection algorithm based on the Abstract Syntax Tree traversal. It helps to explore the uncertainty propagation path (UPP), and assists to transform the deterministic code (DC) into uncertainty intrinsic code (UIC). Based on the open source 2D implicit heat equation \cite{} solver and 3D explicit heat equation sovler \cite{}, the corresponding UIC can be generated. These codes are used to investigate how to gain performance improvements by using UIC.
