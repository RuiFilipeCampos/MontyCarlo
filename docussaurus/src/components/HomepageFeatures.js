import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Easy to Use',
    Svg: require('../../static/img/python.svg').default,
    description: (
      <>
        MontyCarlo was designed from the ground up to be easily installed and
        used to get your simulation up and running quickly.
      </>
    ),
  }, 
  {
    title: 'Focus on What Matters',
    Svg: require('../../static/img/python.svg').default,
    description: (
      <>
        MontyCarlo lets you focus on your simulation, and we&apos;ll do the boring physics stuff.
      </>
    ),
  },
  {
    title: 'Powered by Cython',
    Svg: require('../../static/img/cython.svg').default,
    description: (
      <>
        Extend or customize your website layout by reusing React. Docusaurus can
        be extended while reusing the same header and footer.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} alt={title} />
      </div>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
