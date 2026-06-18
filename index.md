---
title: Home
nav_order: 1
---

<style>
.rem-hero { padding: 2rem 0 1.5rem; border-bottom: 1px solid #e5e7eb; margin-bottom: 2rem; }
.rem-hero h2 { font-size: 1.6rem; font-weight: 700; color: #111827; margin: 0 0 0.75rem; }
.rem-hero p { font-size: 1rem; color: #4b5563; line-height: 1.7; max-width: 680px; margin: 0; }
.rem-grid-top { display: grid; grid-template-columns: repeat(3, 1fr); gap: 16px; padding: 1.5rem 0 0; align-items: stretch; }
.rem-grid-bot { display: grid; grid-template-columns: repeat(2, 1fr); gap: 16px; padding: 0 0 1.5rem; width: 66.6%; margin: 16px auto 0; align-items: stretch; }
.rem-box { background: #fff; border: 1px solid #e5e7eb; border-radius: 12px; padding: 1.25rem; text-decoration: none; display: flex; flex-direction: column; gap: 6px; box-sizing: border-box; }
.rem-title { font-size: 13px; font-weight: 700; color: #111827; margin: 0; line-height: 1.4; }
.rem-sub { font-size: 11px; color: #6b7280; margin: 0; line-height: 1.5; }
.rem-box:hover { filter: brightness(0.97); }
.rem-tag { font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; padding: 3px 8px; border-radius: 6px; display: inline-block; width: fit-content; }
.rem-arrow { font-size: 16px; color: #9ca3af; margin-top: auto; }
.tag-yellow { background: #fef9c3; color: #854d0e; }
.tag-blue   { background: #dbeafe; color: #1e3a8a; }
.tag-red    { background: #ffe4e6; color: #9f1239; }
.tag-green  { background: #dcfce7; color: #14532d; }
.tag-gray   { background: #f3f4f6; color: #374151; }
</style>

## Welcome to this Relational (Hyper) Event Model tutorial space!

This site collects teaching materials on **Relational (Hyper) Event Models** (RHEMs) and their extensions. Whether you are new to R(H)EMs -- with or without H -- or looking to deepen your understanding, you can find here **theory, hands-on tutorials, software, and reading material**. 

<div class="rem-grid-top">

<a class="rem-box" href="00-Notes&Slides/">
<span style="font-size:22px">🗂️</span>
<span class="rem-tag tag-yellow">Notes</span>
<span class="rem-title">Formulary, Slides, and Notes</span>
<span class="rem-sub">Supporting materials for the workshops</span>
<span class="rem-arrow">→</span>
</a>

<a class="rem-box" href="01-Introduction-to-REM/">
<span style="font-size:22px">📈</span>
<span class="rem-tag tag-blue">Workshop 1 - Practical</span>
<span class="rem-title">Introduction to Relational Event Models</span>
<span class="rem-sub">Simulating and fitting REMs</span>
<span class="rem-arrow">→</span>
</a>

<a class="rem-box" href="02-All-you-need-to-know-RHEM/">
<span style="font-size:22px">📊</span>
<span class="rem-tag tag-red">Workshop 2 - Practical</span>
<span class="rem-title">All You Need to Know About Relational Hyper-Event Modeling</span>
<span class="rem-sub">Fitting and interpreting RHEMs</span>
<span class="rem-arrow">→</span>
</a>

</div>

<div class="rem-grid-bot">

<a class="rem-box" href="https://franciscorichter.github.io/amorem/">
<span style="font-size:22px">💻</span>
<span class="rem-tag tag-green">Software</span>
<span class="rem-title"><code>amore</code> package</span>
<span class="rem-sub">R package</span>
<span class="rem-arrow">→</span>
</a>

<a class="rem-box" href="https://arxiv.org/abs/2604.07063">
<span style="font-size:22px">📖</span>
<span class="rem-tag tag-gray">Resource</span>
<span class="rem-title">Introduction to Relational Event Modelling</span>
<span class="rem-sub">Tutorial paper</span>
<span class="rem-arrow">→</span>
</a>

</div>

<div style="margin-top: 2rem; padding: 1.5rem 0; border-top: 1px solid #e5e7eb;">
<p style="font-size: 0.95rem; color: #6b7280; margin: 0;">These materials are a work in progress and we welcome any feedback, suggestions, or just a hello -- feel free to reach out at <a href="mailto:martina.boschi@usi.ch" style="color: #1e3a8a;">martina.boschi@usi.ch</a> or <a href="mailto:melania.lembo@usi.ch" style="color: #1e3a8a;">melania.lembo@usi.ch</a>!</p>
</div>

<!-- ---
title: Home
nav_order: 1
---

## Welcome!

In this site, you can find two tutorials on Relational (Hyper) Event Modeling.

### [Part 1 — Introduction to Relational Event Models](01-Introduction-to-REM/)

Many sociological processes unfold as sequences of temporally ordered interaction events between entities. Leveraging the availability of time-stamped data, the Relational Event Model (REM) provides a powerful framework for generative modeling of such dynamic networks. REMs capture how node and edge characteristics, together with past interactions, shape the evolution of these processes. The essential temporal information embedded in relational event sequences has motivated advances in the specification of covariates and development of inference techniques.

In this workshop, we will introduce the foundations of REMs and show how to construct covariates that represent exogenous drivers, endogenous mechanisms, and temporal features of these processes. We will then examine likelihood-based inference methods for estimating covariate effects and discuss extensions that enhance the flexibility of REMs, including non-linear and time-varying influences as well as random effects to account for network and actor heterogeneity.

### [Part 2 — All You Need to Know About Relational Hyper-Event Modeling](02-All-you-need-to-know-RHEM/)

Social network analysis has relied mainly on time-aggregated data and individual node connections. However, advances in data collection technologies now provide large amounts of time-stamped and hypernetwork information. Examples of such datasets include email multicast communication data and citation networks. Whereas methodological developments in temporal network science have expanded the toolkit for describing time-stamped and hypergraph data, relational hyper-event modeling has emerged as a central approach for the statistical modeling of the dynamic structures these data reveal. Recent contributions have made the Relational Hyper Event Model (RHEM) even more flexible, allowing it to address a range of applied demands.

In this workshop, we aim to present RHEMs, beginning with their core properties and progressing to their more novel modeling capacities, highlighting how they can be employed to address empirical questions. -->