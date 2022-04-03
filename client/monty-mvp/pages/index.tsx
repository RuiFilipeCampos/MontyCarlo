import type { NextPage } from "next";
import { NextRouter, useRouter } from "next/router";

import * as ch from "@chakra-ui/react";

const Home: NextPage = () => {
  const router: NextRouter = useRouter();

  return (
    <>
      <ch.Flex w="100vw" h="100vh">
        <ch.Flex shadow="md" w="full" h="90%" p="10" m="50">
          <ch.VStack>
            <ch.Heading>Monty Carlo Live Demo !</ch.Heading>

            <ch.Heading>Intro</ch.Heading>

            <ch.Text>
              Monty Carlo is a general purpose particle simulator.
            </ch.Text>

            <ch.Text>
              General purpose means that you can construct any geometry and fill
              it with any material !
            </ch.Text>

            <ch.Text>
              Monty already has a lot of features, but it is still very much in
              development.
            </ch.Text>

            <ch.Text>
              This demo is to showcase a fraction of what this program in
              capable, and, for that matter, what codes that are in the category
              "general purpose particle simulators" are capable.
            </ch.Text>

            <ch.Text>
              Other examples are: PENELOPE(the golden standard of Medical
              Physics), GEANT4 (used at CERN), MCNPX (the OG, developed for the
              Manhaten project), [MORE]
            </ch.Text>

            <ch.Button
              variant="solid"
              colorScheme="facebook"
              onClick={() => router.push("/demo")}
            >
              Try Monty Carlo !
            </ch.Button>
          </ch.VStack>
        </ch.Flex>
      </ch.Flex>
    </>
  );
};

export default Home;
