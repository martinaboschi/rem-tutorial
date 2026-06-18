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

---
title: Home
nav_order: 1
---

<style>
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 16px; padding: 1.5rem 0; }
.box { background: #fff; border: 1px solid #e5e7eb; border-radius: 12px; padding: 1.5rem 1.25rem; text-decoration: none; display: flex; flex-direction: column; gap: 10px; transition: border-color 0.15s, background 0.15s; }
.box:hover { border-color: #9ca3af; background: #f9fafb; }
.box-icon { font-size: 24px; color: #6b7280; }
.box-tag { font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; padding: 3px 8px; border-radius: 6px; display: inline-block; width: fit-content; }
.tag-purple { background: #EEEDFE; color: #3C3489; }
.tag-teal   { background: #E1F5EE; color: #085041; }
.tag-coral  { background: #FAECE7; color: #712B13; }
.tag-blue   { background: #E6F1FB; color: #0C447C; }
.box-title { font-size: 15px; font-weight: 600; color: #111827; margin: 0; line-height: 1.4; }
.box-sub   { font-size: 13px; color: #6b7280; margin: 0; line-height: 1.5; }
.arrow { font-size: 18px; color: #9ca3af; margin-top: auto; padding-top: 8px; }
</style>

<div class="grid">

  <a class="box" href="#">
    <span class="box-icon">📖</span>
    <span class="box-tag tag-purple">Workshop 1</span>
    <p class="box-title">Theory</p>
    <p class="box-sub">Foundations of Relational Event Models</p>
    <span class="arrow">→</span>
  </a>

  <a class="box" href="#">
    <span class="box-icon">📖</span>
    <span class="box-tag tag-teal">Workshop 2</span>
    <p class="box-title">Theory</p>
    <p class="box-sub">Foundations of Relational Hyper-Event Models</p>
    <span class="arrow">→</span>
  </a>

  <a class="box" href="01-Introduction-to-REM/">
    <span class="box-icon">💻</span>
    <span class="box-tag tag-purple">Workshop 1</span>
    <p class="box-title">Introduction to Relational Event Models</p>
    <p class="box-sub">Covariates, likelihood inference &amp; extensions</p>
    <span class="arrow">→</span>
  </a>

  <a class="box" href="02-All-you-need-to-know-RHEM/">
    <span class="box-icon">💻</span>
    <span class="box-tag tag-teal">Workshop 2</span>
    <p class="box-title">All You Need to Know About Relational Hyper-Event Modeling</p>
    <p class="box-sub">Core properties &amp; novel modeling capacities</p>
    <span class="arrow">→</span>
  </a>

  <a class="box" href="#">
    <span class="box-icon">📦</span>
    <span class="box-tag tag-coral">Resources</span>
    <p class="box-title">Package</p>
    <p class="box-sub">Software &amp; implementation tools</p>
    <span class="arrow">→</span>
  </a>

  <a class="box" href="#">
    <span class="box-icon">📄</span>
    <span class="box-tag tag-blue">Resources</span>
    <p class="box-title">Tutorial Paper</p>
    <p class="box-sub">Step-by-step methodological guide</p>
    <span class="arrow">→</span>
  </a>

</div>