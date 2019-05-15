## SeuratPlotly
[Plot_ly](https://plot.ly/)-based functions that are enhanced counterparts to the plotting functions available in the [Seurat package](https://github.com/satijalab/seurat).  The primary advantage SeuratPlotly has over the standard plotting functions of Seurat are the inclusion of 3D scatterplots of dimentional reductions.  For example, DimPlotly3D allows viewing the first 3 UMAP dimensions of of the 
[Villani dataset](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5775029/)
<div>
    <a href="https://plot.ly/~milescsmith/1/?share_key=zKLzhDBe2mTLL2bXpSF5cF" target="_blank" title="villani_umap_dimplotly3d" style="display: block; text-align: center;"><img src="https://plot.ly/~milescsmith/1.png?share_key=zKLzhDBe2mTLL2bXpSF5cF" alt="villani_umap_dimplotly3d" style="max-width: 100%;width: 750px;"  width="750" onerror="this.onerror=null;this.src='https://plot.ly/404.png';" /></a>
    <script data-plotly="milescsmith:1" sharekey-plotly="zKLzhDBe2mTLL2bXpSF5cF" src="https://plot.ly/embed.js" async></script>
</div> (click to interact with the plot) whereas FeaturePlotly3D and Feature2Plotly3D allow for viewing the expression of given features (in this case CD14 and CD16) in 3 dimensions:
<div>
    <a href="https://plot.ly/~milescsmith/3/?share_key=HV3YP8852vFIInj7eltyae" target="_blank" title="villani_umap_feature2plotly3d" style="display: block; text-align: center;"><img src="https://plot.ly/~milescsmith/3.png?share_key=HV3YP8852vFIInj7eltyae" alt="villani_umap_feature2plotly3d" style="max-width: 100%;width: 750px;"  width="750" onerror="this.onerror=null;this.src='https://plot.ly/404.png';" /></a>
    <script data-plotly="milescsmith:3" sharekey-plotly="HV3YP8852vFIInj7eltyae" src="https://plot.ly/embed.js" async></script>
</div>
